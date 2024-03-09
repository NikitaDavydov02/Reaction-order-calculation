using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Reaction_orders
{
    public class ReactionSystem
    {
        //Constants
        double k = 1.38 * Math.Pow(10, -23);//Boltzman constant
        double K = Math.Pow(10, -10);
        double E0 = 400;
        double h = 6.62 * Math.Pow(10, -34);
        double Nav = 6.02 * Math.Pow(10, 23);
        double R = 8.314;
        double pi = Math.PI;
        double z = 6;

        //Variables
        public double V { get; private set; }
        public double T { get; private set; }
        public double P { get; private set; }
        public double[] n { get; private set; }
        public double[] gamma { get; private set; }

        double[] chemPot;
        double[] standartChemPot;
        double[] M;
        double[] Q;
        double[][] LJparameters;
        double[][] sigmaMatrix;
        double[][] epsilonMatrix;
        double[][] w;
        double[][] correlator;
        double[][] noncorelatedPairsNumber;
        double[][] corelatedPairsNumber;
        int numberOfComponents = 3;
        double[] m;
        double[] N;
        double[] molarFractions;
        double a;
        double Nsum;

        private ReactionSystem()
        {
            Console.WriteLine("Initializing...");
            M = new double[numberOfComponents];
            n = new double[numberOfComponents];
            Q = new double[numberOfComponents];
            LJparameters = new double[numberOfComponents][];
            sigmaMatrix = new double[numberOfComponents][];
            epsilonMatrix = new double[numberOfComponents][];
            corelatedPairsNumber = new double[numberOfComponents][];
            noncorelatedPairsNumber = new double[numberOfComponents][];
            ////////////////////////////////////
            m = new double[numberOfComponents];
            N = new double[numberOfComponents];
            gamma = new double[numberOfComponents];
            chemPot = new double[numberOfComponents];
            standartChemPot = new double[numberOfComponents];
            w = new double[numberOfComponents][];
            correlator = new double[numberOfComponents][];
            molarFractions = new double[numberOfComponents];
            for (int i = 0; i < numberOfComponents; i++)
            {
                M[i] = 0.03995;
                gamma[i] = 1;
                Q[i] = 1;
                n[i] = 10000;
                LJparameters[i] = new double[2];

                sigmaMatrix[i] = new double[numberOfComponents];
                epsilonMatrix[i] = new double[numberOfComponents];
                corelatedPairsNumber[i] = new double[numberOfComponents];
                noncorelatedPairsNumber[i] = new double[numberOfComponents];
                w[i] = new double[numberOfComponents];
                correlator[i] = new double[numberOfComponents];
            }
            //Argon
            //LJ parameters sigma
            LJparameters[0][1] = 3.40984E-10;
            LJparameters[1][1] = 3.40984E-10;
            LJparameters[2][1] = 3.40984E-10;
            //LJ parameters epsilon
            LJparameters[0][0] = 1.65517E-21;
            LJparameters[1][0] = 1.65517E-21;
            LJparameters[2][0] = 1.65517E-21;
            //Ideal
            for (int i = 0; i < numberOfComponents; i++)
                for (int j = 0; j < numberOfComponents; j++)
                {
                    sigmaMatrix[i][j] = (LJparameters[i][1] + LJparameters[j][1]) / 2.0;
                    epsilonMatrix[i][j] = Math.Sqrt(LJparameters[i][0] * LJparameters[j][0]);
                    w[i][j] = LJpotential(a, i, j) - (LJpotential(a, i, i)+ LJpotential(a, j, j)) / 2;
                }

           // RecalculateCorrelations();
        }
        public ReactionSystem(double _V, double _T, double[] stableReagentsAmounts) : this()
        {
            T = _T;
            for (int i = 0; i < numberOfComponents; i++)
                n[i] = stableReagentsAmounts[i];
            n[2] = K * n[0] * n[1];
            SetVolume(_V);
            SetComponentAmounts(n);
            ReachEquilibriumForTransitionalComplex();
            SetVolume(_V);
        }
        public ReactionSystem(double _P, double[] stableReagentsAmounts, double _T) : this()
        {
            T = _T;
            for (int i = 0; i < numberOfComponents; i++)
                n[i] = stableReagentsAmounts[i];
            if(n[0]!=0 && n[1]!=0)
                n[2] = K * n[0] * n[1];
            SetVolume(2);
            SetComponentAmounts(n);
            ReachEquilibriumForTransitionalComplex();
            SetPressure(_P,FreedomDegree.Volume);

        }
        public void SetComponentAmounts(double[] newAmountsMol)
        {
            Nsum = 0;
            for (int i = 0; i < numberOfComponents; i++)
            {
                N[i] = n[i] * Nav;
                Nsum += N[i];
                m[i] = M[i] / Nav;
            }
            for (int i = 0; i < numberOfComponents; i++)
            {
                molarFractions[i] = N[i] / Nsum;
            }
            
            RecalculateCorrelations();

            /*P = CalculatePressure();
            CalculateStandartChemicalPotentials();
            CalculateChemicalPotentials();
            CalculateActivityCoefficients();*/
        }
        public void SetPressure(double pressure, FreedomDegree freedomDegree)
        {
            if (freedomDegree == FreedomDegree.Volume)
            {
                double targetV = FindVolumeCorrespondingToParticularPressure(pressure);
                SetVolume(targetV);
            }
            else if (freedomDegree == FreedomDegree.N0)
            {
                double target_n0 = FindAmountCorrespondingToParticularPressure(pressure, 0);
                n[0] = target_n0;
                SetComponentAmounts(n);
            }
            else if (freedomDegree == FreedomDegree.N1)
            {
                double target_n1 = FindAmountCorrespondingToParticularPressure(pressure, 1);
                n[1] = target_n1;
                SetComponentAmounts(n);
            }
        }
        public void SetVolume(double volume)
        {
            V = volume;
            a = Math.Pow(3.0 * V / (4.0 * pi * Nsum), 0.33);

            for (int i = 0; i < numberOfComponents; i++)
                for (int j = 0; j < numberOfComponents; j++)
                    w[i][j] = LJpotential(a, i, j) - (LJpotential(a, i, i) + LJpotential(a, j, j)) / 2;

            RecalculateCorrelations();

            /*P = CalculatePressure();
            CalculateStandartChemicalPotentials();
            CalculateChemicalPotentials();
            CalculateActivityCoefficients();*/
        }
        public double CalculateCurrentReactionRate()
        {
            CalculateActivityCoefficients();
            return (k * T / h) * K * gamma[0] * gamma[1] * (N[0] / V) * (N[1] / V) / gamma[2];
        }
        public double GetPressure()
        {
            double dV = 0.01 * V;
            double initV = V;
            SetVolume(initV - 2 * dV);
            double F0 = -k * T * Calculate_LnZ();
            SetVolume(initV - dV);
            double F1 = -k * T * Calculate_LnZ();
            SetVolume(initV + dV);
            double F2 = -k * T * Calculate_LnZ();
            SetVolume(initV + 2 * dV);
            double F3 = -k * T * Calculate_LnZ();
            P = -(8 * F2 - 8 * F1 + F0 - F3) / (12 * dV);
           // P = (F2-F1) / dV;
            SetVolume(initV);
            return P;
        }
        public double Get_F()
        {
            return (-k * T * Calculate_LnZ());
        }
        //--------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------//
        //------------------------------------PRIVATE-------------------------------------------//
        //--------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------//
        private void RecalculateCorrelations()
        {
            double totalNumberOfPairs = z * Nsum / 2;
            double checkSum = 0;
            for (int i = 0; i < numberOfComponents; i++)
                for (int j = 0; j < numberOfComponents; j++)
                {
                    correlator[i][j] = 1;
                    if (i != j)
                        noncorelatedPairsNumber[i][j] = z * N[i] * N[j] / Nsum;
                    else
                        noncorelatedPairsNumber[i][j] = z * N[i] * N[j] / (2 * Nsum);
                    if (j >= i)
                        checkSum += noncorelatedPairsNumber[i][j];
                    corelatedPairsNumber[i][j] = noncorelatedPairsNumber[i][j] * correlator[i][j];
                }
        }
        private void CalculateActivityCoefficients()
        {
            CalculateStandartChemicalPotentials();
            CalculateChemicalPotentials();
            for(int i = 0; i < numberOfComponents; i++)
            {
                if (n[i] == 0)
                {
                    gamma[i] = 1;
                    continue;
                }
                gamma[i] = Math.Exp((chemPot[i] - standartChemPot[i])/(R*T)) / molarFractions[i];
            }
        }
        private void ReachEquilibriumForTransitionalComplex()
        {
            if (n[0] == 0 || n[1] == 0)
                return;
            n[2] = n[0] * n[1] * K * (gamma[0] * gamma[1] / gamma[2]);
            CalculateActivityCoefficients();
            n[2] = n[0] * n[1] * K * (gamma[0] * gamma[1] / gamma[2]);
        }
        private double LnG()
        {
            double output = MathFunctions.StirlingFormula(Nsum);
            for (int i = 0; i < numberOfComponents; i++)
            {
                if (N[i] != 0)
                    output -= MathFunctions.StirlingFormula(N[i]);
            }
            //Noncorrelated pair number
            for (int i = 0; i < numberOfComponents; i++)
            {
                output += MathFunctions.StirlingFormula(noncorelatedPairsNumber[i][i]);
                for (int j = 0; j < numberOfComponents; j++)
                    output += MathFunctions.StirlingFormula(noncorelatedPairsNumber[i][j] / 2);
            }
            for (int i = 0; i < numberOfComponents; i++)
            {
                output -= MathFunctions.StirlingFormula(corelatedPairsNumber[i][i]);
                for (int j = 0; j < numberOfComponents; j++)
                    output -= MathFunctions.StirlingFormula(corelatedPairsNumber[i][j] / 2);
            }
            return output;
        }
        private double Calculate_LnZ()
        {
            // Console.WriteLine("Calculating lnZ");
            double output = 0;
            //log_writer.WriteLine("Calculating lnZ");
            for (int i = 0; i < numberOfComponents; i++)
            {
                if (n[i] == 0)
                    continue;
                //log_writer.WriteLine("    //--------------Component------------");
                double trans = CalculateLnTransitionalStatisticalSum(i);
                //isotherm_writer_line += trans + ";";
                output += trans;
                //log_writer.WriteLine("    LnQtrans = " + trans);
                double freeCellVolume = CalculateFreeCellVolume(i);
                ///isotherm_writer_line += freeCellVolume + ";";
                output += N[i] * Math.Log(freeCellVolume);
                //log_writer.WriteLine("    Free cell volume (m3) = " + freeCellVolume);

                double selfEnergy = z * N[i] * LJpotential(a, i, i) / 2.0;
               // log_writer.WriteLine("    Self energy (J) = " + selfEnergy);
                //log_writer.WriteLine("    Self energy/kT = " + selfEnergy / (k * T));
                //isotherm_writer_line += (-selfEnergy / (k * T)) + ";";
                output -= selfEnergy / (k * T);
            }
            //if (idealGas)
            //    return output;
            for (int i = 0; i < numberOfComponents; i++)
                for (int j = 0; j < i; j++)
                {
                    if (i == j || n[i] == 0 || n[j] == 0)
                        continue;
                    //log_writer.WriteLine("    Interaction energy between " + i + " and " + j);
                    double interactionEnergy = w[i][j] * corelatedPairsNumber[i][j];
                   // log_writer.WriteLine("    Interaction energy (J) =  " + interactionEnergy);
                   // log_writer.WriteLine("    Interaction energy/kT =  " + interactionEnergy / (k * T));
                    output -= interactionEnergy / (k * T);
                   // isotherm_writer_line += (-interactionEnergy / (k * T)) + ";";
                }
            double Ln_G = LnG();
            //isotherm_writer_line += Ln_G + ";";
            output += Ln_G;
           // log_writer.WriteLine("LnG =  " + Ln_G);
            //Calculate free cell volumes
            //Calculate interaction energies
            return output;
        }

        private void CalculateStandartChemicalPotentials()
        {
            double pressure = GetPressure();
            for (int i = 0; i < numberOfComponents; i++)
            {
                double[] cleanCompositions=new double[numberOfComponents];
                for (int j = 0; j < numberOfComponents; j++)
                    cleanCompositions[j] = 0;
                cleanCompositions[i] = 35000;
                double dn = cleanCompositions[i] * 0.01;
                double initial_n = cleanCompositions[i];

                cleanCompositions[i] = initial_n - 2 * dn;
                ReactionSystem system = new ReactionSystem(P,cleanCompositions,T);
                double F0 = system.Get_F();

                cleanCompositions[i] = initial_n - dn;
                system = new ReactionSystem(P, cleanCompositions, T);
                double F1 = system.Get_F();

                cleanCompositions[i] = initial_n + dn;
                system = new ReactionSystem(P, cleanCompositions, T);
                double F2 = system.Get_F();

                cleanCompositions[i] = initial_n + 2 * dn;
                system = new ReactionSystem(P, cleanCompositions, T);
                double F3 = system.Get_F();

                standartChemPot[i] = (8*F2-8*F1+F0-F3) /(12*dn);
            }
        }
        private void CalculateChemicalPotentials()
        {
            double pressure = GetPressure();
            for (int i = 0; i < numberOfComponents; i++)
            {
                if (n[i] == 0)
                {
                    chemPot[i] = 0;
                    continue;
                }
                double[] compositions = new double[numberOfComponents];
                for (int j = 0; j < numberOfComponents; j++)
                    compositions[j] = n[j];
                double dn = compositions[i] * 0.01;
                double initial_n = compositions[i];

                compositions[i] = initial_n - 2 * dn;
                ReactionSystem system = new ReactionSystem(pressure, compositions, T);
                double F0 = system.Get_F();

                compositions[i] = initial_n - dn;
                system = new ReactionSystem(pressure, compositions, T);
                double F1 = system.Get_F();

                compositions[i] = initial_n + dn;
                system = new ReactionSystem(pressure, compositions, T);
                double F2 = system.Get_F();

                compositions[i] = initial_n + 2 * dn;
                system = new ReactionSystem(pressure, compositions, T);
                double F3 = system.Get_F();

                chemPot[i] = (8 * F2 - 8 * F1 + F0 - F3) / (12 * dn);
            }
        }
        private double FindVolumeCorrespondingToParticularPressure(double targetPressure, double Vaccuracy = 0.0000001)
        {
            double Vmin = 0.1;
            double Vmax = 50;

            return MathFunctions.BisectionSolve(PressureEquation, Vaccuracy, Vmin, Vmax, 1, new List<double>() { targetPressure, n[0], n[1], n[2] });

        }
        private double FindAmountCorrespondingToParticularPressure(double targetPressure, int component,double n_accuracy = 0.0000001)
        {
            double min = 0.1*n[component];
            double max= n[component] + 0.1 * n[component];
            return MathFunctions.BisectionSolve(PressureEquation, n_accuracy, min, max, component+2, new List<double>() { targetPressure, n[0], n[1], n[2] });
        }
        private double PressureEquation(List<double> args)
        {
            double targetPressure = args[0];
            double newV = args[1];
            double[] new_n = new double[3];
            new_n[0] = args[2];
            new_n[1] = args[3];
            new_n[2] = args[4];
            SetVolume(newV);
            SetComponentAmounts(new_n);
            GetPressure();

            double output = P - targetPressure;

            return output;
        }
        private double XiFunction(double r, double component)
        {
            double integral = MathFunctions.DefiniteIntegral(Xi_SubintegrativeExpression, -a, a, 0, new List<double>() { r, component });
            double output = integral / (2.0 * a);
            int y = 0;
            return output;
        }
        private double Xi_SubintegrativeExpression(List<double> args)
        {
            //return 0;
            double output = 0;
            double x = args[0];
            double r = args[1];
            int component = (int)args[2];
            for (int i = 0; i < numberOfComponents; i++)
            {
                output += molarFractions[i] * LJpotential(Math.Sqrt(a * a - 2.0 * x * r + r * r), component, i);
            }
            return z * output;
        }
        private double LJpotential(double r, int component1, int component2)
        {
            return 0;
            double sigma_r = sigmaMatrix[component1][component2] / r;

            double output = 4.0 * epsilonMatrix[component1][component2] * (Math.Pow(sigma_r, 12.0) - Math.Pow(sigma_r, 6.0));
            return output;
        }
        private double CalculateFreeCellVolume_SubintegrativeExpression(List<double> args)
        {
            double r = args[1];
            double component = args[0];
            double xi = XiFunction(r, component);
            double output = r * r * Math.Exp(-xi / (k * T));
            if (double.IsNaN(output))
                output = output;
            return output;
        }
        private double CalculateFreeCellVolume(int component)
        {
            if (N[component] == 0)
                throw new Exception();
            double xi = XiFunction(0.0, component);
            double output = 4.0 * pi * Math.Exp(xi / (k * T));

            double volume = MathFunctions.DefiniteIntegral(CalculateFreeCellVolume_SubintegrativeExpression, 0, 0.99 * a, 1, new List<double>() { component });
            if (double.IsNaN(volume))
                output = output;
            output *= volume;
            return output;
        }
        private double CalculateLnTransitionalStatisticalSum(int component)
        {
            if (N[component] == 0)
                throw new Exception();
            return (3.0 * N[component] / 2.0) * Math.Log(2.0 * pi * m[component] * k * T / (h * h));
        }
       
    }
    public enum FreedomDegree
    {
        Volume,
        N0,
        N1,

    }
}
