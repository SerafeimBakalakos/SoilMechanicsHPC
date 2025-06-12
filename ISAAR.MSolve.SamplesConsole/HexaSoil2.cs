using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ISAAR.MSolve.FEM.Entities.TemporalFunctions;

namespace ISAAR.MSolve.SamplesConsole
{
    public class HexaSoil2
    {
        public static double startX = 0.0;
        public static double startY = 0.0;
        public static double startZ = 0.0;
        public static double LengthX = 1.0;
        public static double LengthY = 1.0;
        public static double LengthZ = 1.0;
        int nodeID = 1;
        public static double hx = 5.0;
        public static double hy = 5.0;
        public static double hz = 4.0;
        public static int imax = (int)Math.Truncate(hx / LengthX) + 1;
        public static int jmax = (int)Math.Truncate(hy / LengthY) + 1;
        public static int kmax = (int)Math.Truncate(hz / LengthZ) + 1;
        public static double qfail = 0.0;
        public static double ufail = 0.0;
        public static void MakeHexaSoil(Model model, double[] Stoch1, double[] Stoch2, double[] Stoch3, double[] omega, double lambda)
        {
            // xreiazetai na rythmizei kaneis ta megethi me auto to configuration
            int nodeID = 1;
            for (int l = 0; l < kmax; l++)
            {
                for (int k = 0; k < jmax; k++)
                {
                    for (int j = 0; j < imax; j++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node(nodeID, startX + j * LengthX, startY + k * LengthY, startZ + l * LengthZ));
                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < jmax; j++)
            {
                for (int k = 0; k < imax; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                    model.NodesDictionary[nodeID].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
                    nodeID++;
                }
            }
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };
            Element e1;
            int IDhelp = 1;
            for (int ii = 1; ii < model.NodesDictionary.Count + 1; ii++)
            {
                var nodecheck = model.NodesDictionary[ii];
                if (nodecheck.X != hx && nodecheck.Y != hy && nodecheck.Z != hz)
                {
                    const int gpNo = 8;
                    var initialStresses = new double[6];
                    var element1Nodes = new Node[8];
                    element1Nodes[0] = model.NodesDictionary[ii];
                    element1Nodes[1] = model.NodesDictionary[ii + 1];
                    element1Nodes[2] = model.NodesDictionary[ii + 1 + imax];
                    element1Nodes[3] = model.NodesDictionary[ii + imax];
                    element1Nodes[4] = model.NodesDictionary[ii + jmax * imax];
                    element1Nodes[5] = model.NodesDictionary[ii + jmax * imax + 1];
                    element1Nodes[6] = model.NodesDictionary[ii + jmax * imax + 1 + imax];
                    element1Nodes[7] = model.NodesDictionary[ii + jmax * imax + imax];
                    var nodeCoordinates = new double[8, 3];
                    for (int i = 0; i < 8; i++)
                    {
                        nodeCoordinates[i, 0] = element1Nodes[i].X;
                        nodeCoordinates[i, 1] = element1Nodes[i].Y;
                        nodeCoordinates[i, 2] = element1Nodes[i].Z;
                    }
                    var gaussPointMaterials = new KavvadasClays[8];
                    var gaussPointMaterials2 = new DruckerPragerNLH[8];
                    var young = 2.1e5;
                    var poisson = 0.3;
                    var alpha = 1.0;
                    var ksi = 0.02;
                    var gamma = 10;  //effective stress
                    var Htot = hz;
                    //for (int i = 0; i < gpNo; i++)
                       // gaussPointMaterials[i] = new KavvadasClays(young, poisson, alpha, ksi);
                    var elementType1 = new Hexa8Fixed(gaussPointMaterials);
                    var gaussPoints = elementType1.CalculateGaussMatrices(nodeCoordinates);
                    var elementType2 = new Hexa8u8p(gaussPointMaterials); //this because hexa8u8p has not all the fortran that hexa8 has.
                    var perm = new double[8];
                    for (int i = 0; i < gpNo; i++)
                    {
                        var ActualZeta = 0.0;
                        var ActualXi = 0.0;
                        var ActualPsi = 0.0;
                        double[] Coord = new double[3];
                        var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualXi += help[j] * nodeCoordinates[j, 0];
                            ActualPsi += help[j] * nodeCoordinates[j, 1];
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        //gaussPointMaterials[i].Zeta = ActualZeta;
                        Coord[0] = ActualXi;
                        Coord[1] = ActualPsi;
                        Coord[2] = ActualZeta;
                        var comp = 0.0;
                        ActualZeta = 0.0;
                        for (int j = 0; j < gpNo; j++)
                        {
                            ActualZeta += help[j] * nodeCoordinates[j, 2];
                        }
                        for (int j = 0; j < 8; j = j + 2)
                        {
                            comp += Stoch1[j] * Math.Cos(omega[j] * ActualZeta);
                        }
                        for (int j = 1; j < 8; j = j + 2)
                        {
                            comp += Stoch1[j] * Math.Sin(omega[j] * ActualZeta);
                        }
                        comp = Math.Abs((comp) * 0.25 * 0.008686 + 0.008686) / 1;
                        var csl = 0.0;
                        for (int j = 0; j < 8; j = j + 2)
                        {
                            csl += Stoch2[j] * Math.Cos(omega[j] * ActualZeta);
                        }
                        for (int j = 1; j < 8; j = j + 2)
                        {
                            csl += Stoch2[j] * Math.Sin(omega[j] * ActualZeta);
                        }
                        csl = Math.Abs((csl) * 0.095 * 0.733609251 + 0.733609251) / 1;
                        initialStresses[2] = -gamma * (Htot - ActualZeta);
                        initialStresses[0] = -50;
                        initialStresses[1] = -50;
                        initialStresses[3] = 0;
                        initialStresses[4] = 0;
                        initialStresses[5] = 0;
                        double clayactualzeta = 1.0;
                        if (ActualZeta > clayactualzeta)
                        {
                            //gaussPointMaterials = new KavvadasClays[8];
                            gaussPointMaterials[i] = new KavvadasClays(comp, csl, 1, ksi, initialStresses, Htot, Coord);
                            gaussPointMaterials[i].Zeta = ActualZeta;
                            for (int iiii = 0; iiii < gpNo; iiii++)
                            {
                                var ActualZeta1 = 0.0;
                                var help1 = elementType1.CalcH8Shape(gaussPoints[iiii].Xi, gaussPoints[iiii].Eta, gaussPoints[iiii].Zeta);
                                for (int j = 0; j < gpNo; j++)
                                {
                                    ActualZeta1 += help1[j] * nodeCoordinates[j, 2];
                                }
                                for (int j = 0; j < 8; j = j + 2)
                                {
                                    perm[iiii] += Stoch3[j] * Math.Cos(omega[j] * ActualZeta1);
                                }
                                for (int j = 1; j < 8; j = j + 2)
                                {
                                    perm[iiii] += Stoch3[j] * Math.Sin(omega[j] * ActualZeta1);
                                }
                                perm[iiii] = Math.Abs((perm[iiii]) * 0.25 * Math.Pow(10, -8) + Math.Pow(10, -8)) / 1;
                                perm[iiii] = 3600 * 24 * perm[iiii];
                            }
                            elementType2 = new Hexa8u8p(gaussPointMaterials);
                            //for (int ii = 0; ii < gpNo; ii++)
                            //{
                            //    elementType2.Permeability[ii] = perm[ii];
                            //}
                            for (int j = 0; j < 6; j++)
                            {
                                gaussPointMaterials[i].Stresses[j] = gaussPointMaterials[i].Stresses[j] - gaussPointMaterials[i].initialStresses[j];
                            }
                        }
                        else
                        {
                            //gaussPointMaterials2 = new DruckerPragerNLH[8];
                            double phirad = Math.Asin((3 * Math.Sqrt(1.5) * csl) / (6 + Math.Sqrt(1.5) * csl));
                            gaussPointMaterials2[i] = new DruckerPragerNLH(comp, 1.0 / 3.0, 1, phirad, phirad, "Outer Cone", initialStresses, Coord);
                            for (int iiii = 0; iiii < gpNo; iiii++)
                            {
                                var ActualZeta1 = 0.0;
                                var help1 = elementType1.CalcH8Shape(gaussPoints[iiii].Xi, gaussPoints[iiii].Eta, gaussPoints[iiii].Zeta);
                                for (int j = 0; j < gpNo; j++)
                                {
                                    ActualZeta1 += help1[j] * nodeCoordinates[j, 2];
                                }
                                for (int j = 0; j < 8; j = j + 2)
                                {
                                    perm[iiii] += Stoch3[j] * Math.Cos(omega[j] * ActualZeta1);
                                }
                                for (int j = 1; j < 8; j = j + 2)
                                {
                                    perm[iiii] += Stoch3[j] * Math.Sin(omega[j] * ActualZeta1);
                                }
                                perm[iiii] = Math.Abs((perm[iiii]) * 0.25 * Math.Pow(10, -8) + Math.Pow(10, -8)) / 1;
                                perm[iiii] = 3600 * 24 * perm[iiii] * 1000;
                            }                           
                            elementType2 = new Hexa8u8p(gaussPointMaterials2);
                            //for (int ii = 0; ii < gpNo; ii++)
                            //{
                            //    elementType2.Permeability[ii] = perm[ii];
                            //}
                            for (int j = 0; j < 6; j++)
                            {
                                gaussPointMaterials2[i].Stresses[j] = gaussPointMaterials2[i].Stresses[j] - gaussPointMaterials2[i].initialStresses[j];
                            }
                        }
                    }
                    for (int iii = 0; iii < gpNo; iii++)
                    {
                        elementType2.Permeability[iii] = perm[iii];
                    }
                    //elementType1 = new Hexa8(gaussPointMaterials);
                    //elementType2 = new Hexa8u8p(gaussPointMaterials);
                    //for (int i = 0; i < gpNo; i++)
                    //{
                    //    var ActualZeta = 0.0;
                    //    var help = elementType1.CalcH8Shape(gaussPoints[i].Xi, gaussPoints[i].Eta, gaussPoints[i].Zeta);
                    //    for (int j = 0; j < gpNo; j++)
                    //    {
                    //        ActualZeta += help[j] * nodeCoordinates[j, 2];
                    //    }
                    //    for (int j = 0; j < 8; j = j + 2)
                    //    {
                    //        elementType2.Permeability[i] += Stoch3[j] * Math.Cos(omega[j] * ActualZeta);
                    //    }
                    //    for (int j = 1; j < 8; j = j + 2)
                    //    {
                    //        elementType2.Permeability[i] += Stoch3[j] * Math.Sin(omega[j] * ActualZeta);
                    //    }
                    //    elementType2.Permeability[i] = Math.Abs((elementType2.Permeability[i]) * 0.25 * Math.Pow(10, -8) + Math.Pow(10, -8)) / 1;
                    //}
                    e1 = new Element()
                    {
                        ID = IDhelp
                    };
                    IDhelp++;
                    e1.ElementType = elementType2; //yliko meta diorthosi to material1
                    e1.NodesDictionary.Add(1, model.NodesDictionary[ii]);
                    e1.NodesDictionary.Add(2, model.NodesDictionary[ii + 1]);
                    e1.NodesDictionary.Add(4, model.NodesDictionary[ii + 1 + imax]);
                    e1.NodesDictionary.Add(3, model.NodesDictionary[ii + imax]);
                    e1.NodesDictionary.Add(5, model.NodesDictionary[ii + jmax * imax]);
                    e1.NodesDictionary.Add(6, model.NodesDictionary[ii + jmax * imax + 1]);
                    e1.NodesDictionary.Add(8, model.NodesDictionary[ii + jmax * imax + 1 + imax]);
                    e1.NodesDictionary.Add(7, model.NodesDictionary[ii + jmax * imax + imax]);
                    int subdomainID = 1;
                    //e1.initialForces = e1.ElementType.CalculateForces(e1, initialStresses, initialStresses); no initial forces at this problem we use
                    //for (int i = 0; i < gpNo; i++)
                    //    for (int j = 0; j < 6; j++)
                    //    {
                    //        gaussPointMaterials[i].Stresses[j] = gaussPointMaterials[i].Stresses[j] - gaussPointMaterials[i].initialStresses[j];
                    //    }
                    model.ElementsDictionary.Add(e1.ID, e1);
                    model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
                }
            };
            //DOFType doftype1;
            //doftype1 = new DOFType();
            #region initalloads
            //Load loadinitialx = new Load();
            //Load loadinitialy = new Load();
            //Load loadinitialz = new Load();
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialx = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.X
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //        var help1 = 0;
            //        var help = 0;
            //        foreach(Node nodeelement in elementcheck.NodesDictionary.Values)
            //        {
            //            if (nodeelement==nodecheck)
            //            {
            //                    help = help1;
            //            }
            //                help1 += 1;
            //        }
            //            loadinitialx.Amount += -elementcheck.initialForces[3 * help];
            //        }
            //    }
            //    model.Loads.Add(loadinitialx);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialy = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Y
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            loadinitialy.Amount += -elementcheck.initialForces[3*help+1];
            //        }
            //    }
            //    model.Loads.Add(loadinitialy);
            //}
            //foreach (Node nodecheck in model.NodesDictionary.Values)
            //{
            //    loadinitialz = new Load()
            //    {
            //        Node = nodecheck,
            //        //DOF = doftype1,
            //        DOF = DOFType.Z
            //    };
            //    foreach (Element elementcheck in model.ElementsDictionary.Values)
            //    {
            //        var Pa = -5.0;
            //        var P2a = -10.0/2;
            //        var P4a = -20.0/4;
            //        var bool1 = elementcheck.NodesDictionary.ContainsValue(nodecheck);
            //        var bool2 = nodecheck.Z == hz;
            //        var bool3 = (nodecheck.X == 0 || nodecheck.X == hx) && (nodecheck.Y == 0 || nodecheck.Y == hy);
            //        var bool4 = nodecheck.X != 0 && nodecheck.X != hx && nodecheck.Y != 0 && nodecheck.Y != hy;
            //        if (bool1 == true)
            //        {
            //            var help1 = 0;
            //            var help = 0;
            //            foreach (Node nodeelement in elementcheck.NodesDictionary.Values)
            //            {
            //                if (nodeelement == nodecheck)
            //                {
            //                    help = help1;
            //                }
            //                help1 += 1;
            //            }
            //            if (bool2 == true)
            //            {
            //                if (bool3 == true)
            //                {
            //                    loadinitialz.Amount += -Pa - elementcheck.initialForces[3 * help+2];
            //                }
            //                else if (bool4 == true)
            //                {
            //                    loadinitialz.Amount += -P4a - elementcheck.initialForces[3 * help+2];
            //                }
            //                else
            //                {
            //                    loadinitialz.Amount += -P2a - elementcheck.initialForces[3 * help+2];
            //                }
            //            }
            //            else
            //            {
            //                loadinitialz.Amount += -elementcheck.initialForces[3 * help+2];
            //            }
            //        }
            //    }
            //    model.Loads.Add(loadinitialz);
            //}
            #endregion

            double nodalLoad = 0.0;
            double phideg = 60.0;
            double totalDuration = 1;
            double timeStepDuration = 0.001;
            double constantsegmentdurationratio = 1;
            GeneralDynamicNodalLoad loadinitialz;
            nodalLoad = 0.0;
            foreach (Node nodecheck in model.NodesDictionary.Values)
            {
                nodalLoad = 0.0;
                if (nodecheck.X == 2 && nodecheck.Y == 2 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = -lambda * 1.0 / 4.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Cos(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 2 && nodecheck.Y == 3 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = -lambda * 1.0 / 4.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Cos(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 3 && nodecheck.Y == 2 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = -lambda * 1.0 / 4.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Cos(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 3 && nodecheck.Y == 3 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = -lambda * 1.0 / 4.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Cos(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodalLoad != 0.0)
                {
                    var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, constantsegmentdurationratio * totalDuration);
                    loadinitialz = new GeneralDynamicNodalLoad(nodecheck, StructuralDof.TranslationX, timeFunction);
                    model.TimeDependentNodalLoads.Add(loadinitialz);
                }
            }
            foreach (Node nodecheck in model.NodesDictionary.Values)
            {
                nodalLoad = 0.0;
                if (nodecheck.X == 2 && nodecheck.Y == 2 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = 0.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Sin(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 2 && nodecheck.Y == 3 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = 0.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Sin(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 3 && nodecheck.Y == 2 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = 0.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Sin(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodecheck.X == 3 && nodecheck.Y == 3 && nodecheck.Z == 4)
                {
                    if (phideg == 0.0)
                    {
                        nodalLoad = 0.0;
                    }
                    else
                    {
                        nodalLoad = -lambda * Math.Sin(phideg * Math.PI / 180) / 2.0;
                    }
                }
                if (nodalLoad != 0.0)
                {
                    var timeFunction = new RampTemporalFunction(nodalLoad, totalDuration, timeStepDuration, constantsegmentdurationratio * totalDuration);
                    loadinitialz = new GeneralDynamicNodalLoad(nodecheck, StructuralDof.TranslationZ, timeFunction);
                    model.TimeDependentNodalLoads.Add(loadinitialz);
                }
            }




        }
    }
}


