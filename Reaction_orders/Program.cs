using System.Net.Http.Headers;
using System.Net.NetworkInformation;
using System.IO;
using System.Diagnostics.Tracing;

namespace Reaction_orders
{
    internal class Program
    {
        static double k = 1.38*Math.Pow(10, -23);//Boltzman constant
        static double K = Math.Pow(10, -10);
        static double T =15;//Temperature
        // static double M_a = 18;    //Molar mass
        //static double M_b = 18;    //Molar mass

        // static double r_a =Math.Pow(10, -9) * 0.14;
        //static double r_b =Math.Pow(10, -9)*0.14;
        //static double n_a = 10000;
        //static double n_b = 10000;
        //static double E0 = 40000;//J/mol
        static double E0 = 400;//J/mol
        static double V = 1;//Volume
        //static double Q_a = 1;
        //static double Q_b = 1;
        static double h= 6.62 * Math.Pow(10, -34);
        static double Nav = 6.02 * Math.Pow(10, 23);
        static double R = 8.314;
        static double pi = Math.PI;
        static double z = 6;
        static double[] M;//Molar masses
       // static double[] r;
        static double[] n;//Amounts
        static double[] gamma;
        static double[] Q;
        static double[][] LJparameters;
        static double[][] sigmaMatrix;
        static double[][] epsilonMatrix;
        static double[][] w;
        static double[][] beta;
        static double[][] noncorelatedPairsNumber;
        static double[][] corelatedPairsNumber;
        static int numberOfComponents = 2;
        //Calculate parameters
        static double[] m;
        static double[] N;
        static double[] molarFractions;
        static double a;
        static double Nsum;

        //Calculatin process settings
        static double minV = 2;
        static double maxV = 3;
        //static double maxV = 20;
        static double nVSteps = 100;

        static bool log = false;
        static StreamWriter log_writer;
        static StreamWriter isotherm_writer;
        static StreamWriter chemicalPotential_writer;

        delegate double Function(double x);
        
        static void Main(string[] args)
        {
            log_writer = new StreamWriter(File.Create("log.txt"));
            chemicalPotential_writer = new StreamWriter(File.Create("chemical_potentia.txt"));


            Console.WriteLine("Initializing...");
            M = new double[numberOfComponents];
            //r = new double[numberOfComponents];
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
            w = new double[numberOfComponents][];
            beta = new double[numberOfComponents][];
            molarFractions = new double[numberOfComponents];
            for (int i = 0; i < numberOfComponents; i++)
            {
                M[i] = 0.03995;
                gamma[i] = 1;
               // r[i] = Math.Pow(10, -9) * 0.14;
                Q[i] = 1;
                n[i] = 10000;
                LJparameters[i] = new double[2];
                sigmaMatrix[i] = new double[numberOfComponents];
                epsilonMatrix[i] = new double[numberOfComponents];
                corelatedPairsNumber[i] = new double[numberOfComponents];
                noncorelatedPairsNumber[i] = new double[numberOfComponents];
                w[i] = new double[numberOfComponents];
                beta[i] = new double[numberOfComponents];
            }
            n[0] = 34928.66083 / 2.0;
            n[1] = 34928.66083/2.0;
            //n[2] = 0;
            //Argon
            //epsilon (J/molecule)
            LJparameters[0][0] = 1.65517E-21;
            LJparameters[1][0] = 1.75517E-21;
            //sigma (m)
            LJparameters[0][1] = 3.40984E-10;
            LJparameters[1][1] = 3.30984E-10;
            //Argom end
            //LJ parameters
            /*LJparameters[0][0] = 1.65517E-21;
            LJparameters[1][0] = 1.65517E-21;
            LJparameters[2][0] = 1.65517E-21;

            LJparameters[0][1] = 3.40984E-10;
            LJparameters[1][1] = 3.40984E-10;
            LJparameters[2][1] = 3.40984E-10;*/
            //LJ parameters

            //Ideal gas
            /* //epsilon (J/molecule)
             LJparameters[0][0] = 0;
             LJparameters[1][0] = 0;
             //sigma (m)
             LJparameters[0][1] = 3.40984E-10;
             LJparameters[1][1] = 3.40984E-10;*/


            //Isotherm(true);
            //Isotherm();

            //double bestV = FindVolumeCorrespondingToParticularPressure(101325);
            //Console.WriteLine("Best volume = " + bestV);

            /* double rrate = CalculateCurrentReactionRate() / (Nav*1000); //[mol/(s*l)]
             Console.WriteLine("Reactio rate: " + rrate + " mol/(s*l)");
             Console.WriteLine("Concentration: " + (N[0] / (V * Nav * 1000)) + " M");*/

            //CalculateReactionOrders();

            ReactionOrdersForDifferntCompositions(20);

            //Console.WriteLine("Pressure: " + CalculatePressureAtCurrentConditions());

            log_writer.Close();
            chemicalPotential_writer.Close();
            Console.WriteLine("Main is done");
         }
         static string isotherm_writer_line;
            /// <summary>
            /// Calculates volume for particular pressure and composition of system
            /// </summary>
            /// <param name="targetPressure"></param>
            /// <param name="Vaccuracy"></param>
            /// <returns></returns>
         static double FindVolumeCorrespondingToParticularPressure(double targetPressure, double Vaccuracy=0.001)
         {
            return BisectionSolve(ComponentPartialVolumeEquation, Vaccuracy, 0.1, 50, 1, new List<double>() { targetPressure, n[0], n[1] });
            ////////////////////OLD ALGORITHM//////////////////////////////////////////////////
            double Vstep;
             if (Vaccuracy== -1)
                 Vstep = (maxV - minV) / 100;
             else
                 Vstep = Vaccuracy;
             List<KeyValuePair<double, double>> VvsP = new List<KeyValuePair<double, double>>();
             List<KeyValuePair<double, double>> VvsF = new List<KeyValuePair<double, double>>();

             log_writer.WriteLine("Find Volume Corresponding To ParticularPressure");
             log_writer.WriteLine("targetPressure: "+ targetPressure);
             log_writer.WriteLine("Vmin: " + minV);
             log_writer.WriteLine("Vmax: " + maxV);
             for (V = minV; V <= maxV; V += Vstep)
             {
                 UpdateParameters();
                 double F = -k * T * Calculate_LnZ();
                 log_writer.WriteLine("    V: " + V +"    F: "+ F);
                 VvsF.Add(new KeyValuePair<double, double>(V, F));
                 Console.WriteLine("Volume: " + V + "     F: " + F);
             }
             log_writer.WriteLine("    Pressure calculation");
             for (int i = 1; i < VvsF.Count; i++)
             {
                 double pressure = -(VvsF[i].Value - VvsF[i - 1].Value) / (VvsF[i].Key - VvsF[i - 1].Key);
                 VvsP.Add(new KeyValuePair<double, double>(VvsF[i].Key, pressure));
                 log_writer.WriteLine("    V: " + V + "    P: " + pressure);
                 Console.WriteLine("Volume: " + VvsF[i].Key + "     Pressure: " + pressure);
             }
             int iterationBest = -1;
             for(int i = 0; i < VvsP.Count-1; i++)
             {
                 if (targetPressure - VvsP[i].Value<0 && targetPressure - VvsP[i+1].Value >0)
                 {
                     iterationBest = i;
                     break;
                 }
             }

             if (iterationBest < 0)
             {
                 log_writer.WriteLine("    No best iteration");
                 log_writer.WriteLine("-----------------------------------Find Volume Corresponding To ParticularPressure-------------------------------------");
                 return 0;
             }
             double outputV = VvsP[iterationBest].Key + (VvsP[iterationBest+1].Key - VvsP[iterationBest].Key) * (targetPressure- VvsP[iterationBest].Value) / (VvsP[iterationBest+1].Value - VvsP[iterationBest].Value);
             log_writer.WriteLine("    Best volume: " + outputV);
             log_writer.WriteLine("-----------------------------------Find Volume Corresponding To ParticularPressure-------------------------------------");
             return outputV;
         }
         /// <summary>
         /// Calculates reaction orders for different compositions of two component system
         /// </summary>
         /// <param name="steps"></param>
         /// <param name="n_sum"></param>
        static void ReactionOrdersForDifferntCompositions(int steps = 5, double n_sum=35000)
        {
            StreamWriter orderWriter = new StreamWriter("orders.csv");
            orderWriter.WriteLine("x0;x1;V;absolute;relative");
            double xStep = 1 / (double)steps;
            for(double x0 = 0.01; x0 < 1; x0 += xStep)
            {
                Console.WriteLine("-----------------------NEW COMPOSITION-------------------------");
                n[0] = n_sum * x0;
                n[1] = n_sum * (1 - x0);
               // n[2] = 0;
                V = FindVolumeCorrespondingToParticularPressure(101325);
               // UpdateEquilibriumComposition(n);
              //  chemicalPotential_writer.WriteLine(gamma[0] + "    " + gamma[1] + "    " + gamma[2]);

                
                Console.WriteLine("V=" +V);
                List<double> orders = CalculateReactionOrders();
                Console.WriteLine("Absolute order A: " + orders[0]);
                Console.WriteLine("Relative order A: " + orders[1]);
                Console.WriteLine("Absolute order B: " + orders[2]);
                Console.WriteLine("Relative order B: " + orders[3]);
                orderWriter.WriteLine(x0 + ";" + (1 - x0) + ";" + V + ";" + orders[0] + ";" + orders[1]+ ";" + orders[2] + ";" + orders[3]);
                
            }
            orderWriter.Close();
        }
        static void UpdateEquilibriumComposition(double[] initial_n)
        {
            n[0] = initial_n[0];
            n[1] = initial_n[1];
            n[2] = K * n[0] * n[1]  *gamma[0] * gamma[1] / gamma[2];
            CalculateActivityCoeffitientsForParticularComposition(n);
            n[2] = K * n[0] * n[1] * gamma[0] * gamma[1] / gamma[2];
            n[0] -= n[2];
            n[1] -= n[2];
        }
        static void CalculateActivityCoeffitientsForParticularComposition(double[] composition)
        {
            double[] old_composition = new double[numberOfComponents];
            for (int i = 0; i < numberOfComponents; i++)
            {
                old_composition[i] = n[i];
                n[i] = composition[i];
            }
            double F = -k * T * Calculate_LnZ();    
            for(int i = 0; i < numberOfComponents; i++)
            {
                //Calculating chemical potentials for each component
                double dn = old_composition[i] * 0.01;
                n[i] += dn;
                double F2 = -k * T * Calculate_LnZ();
                n[i] = old_composition[i];
                double _gamma = (F2-F) / dn;
                gamma[i] = _gamma;
            }

        }
        /// <summary>
        /// Calculates all types of reaction orders for particular conditions
        /// </summary>
        /// <returns></returns>
         static List<double> CalculateReactionOrders()
         {
            List<double> output = new List<double>();
            double initialRate = CalculateCurrentReactionRate();
            double initialPressure = CalculatePressureAtCurrentConditions();
            Console.WriteLine("Initial pressure: " + initialPressure);
            Console.WriteLine("Initial reaction rate:" + initialRate);
            double initV = V;
            double init_n0 = n[0];
            double init_n1 = n[1];

            //Absolute reaction order calculation
            n[0] *= 1.01;
            double newPressure = CalculatePressureAtCurrentConditions();
            double delta_n0 = n[0] - init_n0;
            double rateAfterNonIsobaricAddition = CalculateCurrentReactionRate();
            Console.WriteLine("Reaction rate after non-isobaric addition:" + rateAfterNonIsobaricAddition);
            double absoluteOrder = (Math.Log(rateAfterNonIsobaricAddition)-Math.Log(initialRate))/ (Math.Log(n[0]/V)- Math.Log(init_n0 / V));
            Console.WriteLine("Absolute order of A: " + absoluteOrder);

            //Relative order calculation
            //n[0] = init_n0;
            double target_n1 = BisectionSolve(ComponentPartialVolumeEquation, delta_n0 * 0.01, init_n1 - 10 * delta_n0, init_n1 + 10 * delta_n0, 3, new List<double>() { initialPressure, V, n[0] });
            n[1] = target_n1;
            double rateAfterIsobaricAddition = CalculateCurrentReactionRate();
            Console.WriteLine("Reaction rate after isobaric addition: " + rateAfterIsobaricAddition);
            double relativeOrder = (Math.Log(rateAfterIsobaricAddition) - Math.Log(initialRate)) / (Math.Log(n[0] / V) - Math.Log(init_n0 / V));
            Console.WriteLine("Relative order of A: " + relativeOrder);

            n[0] = init_n0;
            n[1] = init_n1;
            V = initV;

            //Absolute reaction order calculation
            n[1] *= 1.01;
            double newPressure1 = CalculatePressureAtCurrentConditions();
            double delta_n1 = n[1] - init_n1;
            double rateAfterNonIsobaricAdditionOfB = CalculateCurrentReactionRate();
            Console.WriteLine("Reaction rate after non-isobaric addition:" + rateAfterNonIsobaricAdditionOfB);
            double absoluteOrderB = (Math.Log(rateAfterNonIsobaricAdditionOfB) - Math.Log(initialRate)) / (Math.Log(n[1] / V) - Math.Log(init_n1 / V));
            Console.WriteLine("Absolute order of B: " + absoluteOrderB);

            //Relative order calculation
            //n[0] = init_n0;
            double target_n0 = BisectionSolve(ComponentPartialVolumeEquation, delta_n1 * 0.01, init_n0 - 10 * delta_n1, init_n0 + 10 * delta_n1, 2, new List<double>() { initialPressure, V, n[1] });
            n[0] = target_n0;
            double rateAfterIsobaricAdditionB = CalculateCurrentReactionRate();
            Console.WriteLine("Reaction rate after isobaric addition of B: " + rateAfterIsobaricAdditionB);
            double relativeOrderofB = (Math.Log(rateAfterIsobaricAdditionB) - Math.Log(initialRate)) / (Math.Log(n[1] / V) - Math.Log(init_n1 / V));
            Console.WriteLine("Relative order of A: " + relativeOrderofB);

            n[0] = init_n0;
            n[1] = init_n1;
            V = initV;

            output.Add(absoluteOrder);
            output.Add(relativeOrder);
            output.Add(absoluteOrderB);
            output.Add(relativeOrderofB);
            return output;
        }
         static double CalculateCurrentReactionRate()
         {
             UpdateParameters();
             double reducedMass = (M[0] * M[1] / (M[0] + M[1]));
             double sqrt = Math.Sqrt(8 * R * T / reducedMass);
             double reactionRate = sqrt * sigmaMatrix[0][1] * sigmaMatrix[0][1] * (N[0]/V) * (N[1]/V) * (2 / (1 + beta[0][1])) * Math.Exp(-E0 / (R * T));
             //[1/(s*m3)]
             return reactionRate;
         }
         static void Isotherm(bool idealGas=false)
         {
             Console.WriteLine("Isotherm calculation start");
             Console.WriteLine("Volume, m^3                   Pressure, Pa");
             double Vstep = (maxV - minV) / nVSteps;
             StreamWriter isotherm_writer = new StreamWriter(File.Create("isotherm.txt"));
             xi_r_writer = new StreamWriter(File.Create("xi_r"));
             isotherm_writer.WriteLine("Volume; a; Vcell;Qtrans1;freeV1; Uself1_kt;Qtrans2;freeV2; Uself2_kt;U12_kT;lnG;F; density");
                 for (V = minV; V <= maxV; V += Vstep)
                 {
                     isotherm_writer_line = V+";";
                     log_writer.WriteLine("-----------------NEW VOLUME-------------------");
                     log_writer.WriteLine("V = " + V);
                     UpdateParameters();
                     double F;
                     if (idealGas)
                         F = -k * T * Calculate_LnZ_Gas();
                     else
                         F = -k * T * Calculate_LnZ();
                     isotherm_writer_line += F + ";";
                     log_writer.WriteLine("F = " + F);
                     double density =( N[0]*M[0]+ N[0] * M[0] )/ (Nav*V);
                     isotherm_writer_line += density + ";";
                     isotherm_writer_line += (LJpotential(a, 0, 0)/(k*T)) + ";";
                     Console.WriteLine("Volume: " + V + " ; F " + F);
                     isotherm_writer.WriteLine(isotherm_writer_line);
                 }
             isotherm_writer.Close();
             xi_r_writer.Close();
             Console.WriteLine("Isotherm calculation is done");
         }
        static double CalculatePressureForParticularVolumeAndAmounts(double newV, double[] amounts)
        {
            V = newV;
            for (int i = 0; i < amounts.Length; i++)
                n[i] = amounts[i];
            double output = CalculatePressureAtCurrentConditions();
            if (double.IsNaN(output))
                V = V;
            return output;

        } 
        static double ComponentPartialVolumeEquation(double[] args)
        {
            double targetPressure = args[0];
            double newV = args[1];
            double new_n0 = args[2];
            double new_n1 = args[3];
            double output = CalculatePressureForParticularVolumeAndAmounts(newV, new double[] { new_n0, new_n1 }) - targetPressure;
            if (double.IsNaN(output))
                output=output;
            return output;
        }
         static double CalculatePressureAtCurrentConditions()
         {
             UpdateParameters();
             double oldF = -k * T * Calculate_LnZ();
             double oldVolume = V;
             //double newVolume = 1.01*V;
             double newVolume = V + 0.01;
             V = newVolume;
             UpdateParameters();
             double newF = -k * T * Calculate_LnZ();
             V = oldVolume;
             UpdateParameters();
             return -(newF-oldF) / (newVolume-oldVolume);
         }
         static double Calculate_LnZ_Gas()
         {
             // Console.WriteLine("Calculating lnZ");
             double output = 0;
             log_writer.WriteLine("Calculating lnZ");
             for (int i = 0; i < numberOfComponents; i++)
             {
                 log_writer.WriteLine("    //--------------Component------------");
                 double trans = CalculateLnTransitionalStatisticalSum(i);
                 isotherm_writer_line += trans + ";";
                 output += trans;
                 log_writer.WriteLine("    LnQtrans = " + trans);
                 output += N[i] * Math.Log(V);

                // double selfEnergy = z * N[i] * LJpotential(a, i, i) / 2.0;
                 //log_writer.WriteLine("    Self energy (J) = " + selfEnergy);
                 //log_writer.WriteLine("    Self energy/kT = " + selfEnergy / (k * T));
                 //isotherm_writer_line += (-selfEnergy / (k * T)) + ";";
                 //output -= selfEnergy / (k * T);
             }
             //if (idealGas)
             //    return output;
             /*for (int i = 0; i < numberOfComponents; i++)
                 for (int j = 0; j < i; j++)
                 {
                     if (i == j)
                         continue;
                     log_writer.WriteLine("    Interaction energy between " + i + " and " + j);
                     double interactionEnergy = (N[i] * N[j]*LJpotential(a,i,j)) / ((N[i] + N[j]));
                     log_writer.WriteLine("    Interaction energy (J) =  " + interactionEnergy);
                     log_writer.WriteLine("    Interaction energy/kT =  " + interactionEnergy / (k * T));
                     output -= interactionEnergy / (k * T);
                     isotherm_writer_line += (-interactionEnergy / (k * T)) + ";";
                 }*/
            return output;
        }
        static double Calculate_LnZ()
        {
           // Console.WriteLine("Calculating lnZ");
            double output = 0;
            log_writer.WriteLine("Calculating lnZ");
            for (int i = 0; i < numberOfComponents; i++)
            {
                log_writer.WriteLine("    //--------------Component------------");
                double trans = CalculateLnTransitionalStatisticalSum(i);
                isotherm_writer_line += trans + ";";
                output += trans;
                log_writer.WriteLine("    LnQtrans = " + trans);
                double freeCellVolume = CalculateFreeCellVolume(i);
                isotherm_writer_line += freeCellVolume + ";";
                output += N[i] * Math.Log(freeCellVolume);
                log_writer.WriteLine("    Free cell volume (m3) = " + freeCellVolume);

                double selfEnergy = z * N[i] * LJpotential(a, i, i) / 2.0;
                log_writer.WriteLine("    Self energy (J) = " + selfEnergy);
                log_writer.WriteLine("    Self energy/kT = " + selfEnergy / (k * T));
                isotherm_writer_line += (-selfEnergy / (k * T)) + ";";
                output -=selfEnergy / (k * T);
            }
            //if (idealGas)
            //    return output;
             for (int i=0;i<numberOfComponents;i++)
                for (int j = 0; j < i; j++)
                {
                    if (i == j)
                        continue;
                    log_writer.WriteLine("    Interaction energy between " + i + " and " + j);
                    double interactionEnergy = (w[i][j] * N[i] * N[j] * 2.0) / ((beta[i][j] + 1.0) * (N[i] + N[j]));
                    log_writer.WriteLine("    Interaction energy (J) =  " + interactionEnergy);
                    log_writer.WriteLine("    Interaction energy/kT =  " + interactionEnergy / (k * T));
                    output -= interactionEnergy / (k * T);
                    isotherm_writer_line += (-interactionEnergy / (k * T)) + ";";
                }
            double Ln_G = LnG();
            isotherm_writer_line +=Ln_G + ";";
            output += Ln_G;
            log_writer.WriteLine("LnG =  " + Ln_G);
            //Calcuate free cell volumes
            //Calculate interaction energies
            return output;
        }
        static double LnG()
        {
            double output = StirlingFormila(Nsum);
            for(int i = 0; i < numberOfComponents; i++)
            {
                output -= StirlingFormila(N[i]);
            }
            output += StirlingFormila(z * (N[0] - noncorelatedPairsNumber[0][1]) / 2.0);
            output += StirlingFormila(z * (N[1] - noncorelatedPairsNumber[0][1]) / 2.0);
            output += StirlingFormila(z * (noncorelatedPairsNumber[0][1]) / 2.0);
            output += StirlingFormila(z * (noncorelatedPairsNumber[0][1]) / 2.0);

            output -= StirlingFormila(z * (N[0] - corelatedPairsNumber[0][1]) / 2.0);
            output -= StirlingFormila(z * (N[1] - corelatedPairsNumber[0][1]) / 2.0);
            output -= StirlingFormila(z * (corelatedPairsNumber[0][1]) / 2.0);
            output -= StirlingFormila(z * (corelatedPairsNumber[0][1]) / 2.0);
            return output;
        }
        static double StirlingFormila(double N)
        {
            return (N * Math.Log(N) - N);
        }
        static double CalculateLnTransitionalStatisticalSum(int component)
        {
           // Console.WriteLine("Calculating tranistional statistical sum");
            return (3.0 * N[component] / 2.0)*Math.Log(2.0 * pi * m[component] * k * T / (h * h));
        }
        static double CalculateFreeCellVolume(int component)
        {
            //if (a > 3.76 * Math.Pow(10, -10))
                log = true;
            double xi = XiFunction(0.0, component);
            double output = 4.0 * pi * Math.Exp(xi / (k * T));
            
            double volume =  DefiniteIntegral(CalculateFreeCellVolume_SubintegrativeExpression, 0, 0.99*a,1,new List<double>() { component});
            if (double.IsNaN(volume))
                output = output;
            output *= volume;
            
            if (log)
                log = false;
            log_writer.WriteLine("-------------------------------------------------");
            log_writer.WriteLine("Calculating free cell volume = "+output);
            log_writer.WriteLine("Cell volume = " + 4.0*pi*a*a*a/3.0);
            log_writer.WriteLine("a = " + a);
            log_writer.WriteLine("-------------------------------------------------");

            //return 4.0 * pi * a * a * a / 3.0;
            return output;
        }
        static double CalculateFreeCellVolume_SubintegrativeExpression(List<double> args)
        {
            double r = args[1];
            double component = args[0];
            double xi =  XiFunction(r, component);
            double output = r * r * Math.Exp(-xi/ (k * T));
            if (double.IsNaN(output))
                output = output;
            return output;
        }
        static void UpdateParameters()
        {
            //Console.WriteLine("Updating parameters");
            Nsum = 0;
            for (int i = 0; i < numberOfComponents; i++)
            {
                N[i] = n[i] * Nav;
                Nsum += N[i];
                m[i] = M[i] / Nav;
            }
            log_writer.WriteLine("///////////////////Updating parameters////////////////////");
            for (int i = 0; i < numberOfComponents; i++)
            {
                molarFractions[i] = N[i] / Nsum;
                log_writer.WriteLine("N "+ i + " = " + N[i]);
            }
                
            a = Math.Pow(3.0 * V / (4.0 * pi * Nsum), 0.33);

            isotherm_writer_line += a + ";" + (4.0*pi*a*a*a/3.0)+";";
            log_writer.WriteLine("a (m) =  " + a);
            

            for (int i=0;i<numberOfComponents;i++)
                for(int j = 0; j < numberOfComponents; j++)
                {
                    sigmaMatrix[i][j] = (LJparameters[i][1] + LJparameters[j][1]) / 2.0;
                    epsilonMatrix[i][j] = Math.Sqrt(LJparameters[i][0] * LJparameters[j][0]);
                    
                }
            for(int i=0;i<numberOfComponents;i++)
                for(int j = 0; j < numberOfComponents; j++)
                {
                    w[i][j] = LJpotential(a, i, j) - (LJpotential(a, i, i) + LJpotential(a, j, j)) / 2.0;
                    beta[i][j] = Math.Sqrt(1 + 4.0 * molarFractions[i] * molarFractions[j] * (Math.Exp(2.0 * w[i][j] / (k * T)) - 1.0));
                    noncorelatedPairsNumber[i][j] = N[i] * N[j] / (N[i] + N[j]);
                    corelatedPairsNumber[i][j] = noncorelatedPairsNumber[i][j] * 2.0 / (1.0 + beta[i][j]);
                }
            OutputMatix("Sigma matrix", sigmaMatrix, numberOfComponents);
            OutputMatix("Epsilon matrix", epsilonMatrix, numberOfComponents);
            OutputMatix("W matrix", w, numberOfComponents);
            OutputMatix("Beta matrix", beta, numberOfComponents);
            GrapXiFromR();
        }
        static void OutputMatix(string name, double[][] matrix, int n)
        {
            log_writer.WriteLine("---------------------------------------");
            log_writer.WriteLine(name);
            for(int i = 0; i < n; i++)
            {
                string line = "";
                for (int j = 0; j < n; j++)
                    line += " " + matrix[i][j];
                log_writer.WriteLine(line);
            }
        }
        static double DefiniteIntegral(Func<List<double>, double> func, double a, double b,int indexOfIntegratingVariable, List<double> nonIntegratingVariablesValues, double nSteps=1000)
        {
            if (indexOfIntegratingVariable == 1)
                ;
            //Console.WriteLine("Itegrating");
            double output = 0;
            double step = (b - a) / (double)nSteps;
            double x_plus_dx;
            double dS;
            List<double> arguments = new List<double>();
            int currentNonIntegrativeVariable = 0;
            for(int i = 0; i <= nonIntegratingVariablesValues.Count; i++)
            {
                if (i == indexOfIntegratingVariable)
                {
                    arguments.Add(a);
                    continue;
                }
                arguments.Add(nonIntegratingVariablesValues[currentNonIntegrativeVariable]);
                currentNonIntegrativeVariable++;
            }
            string indent = "";
            if (log)
            {
                if (a != 0)
                    indent = "      ";
               /* Console.WriteLine(indent+"--------------------Integration----------------------");
                Console.WriteLine(indent + "------------------------------------------------------");
                Console.WriteLine(indent + "------------------------------------------------------");
                Console.WriteLine(indent + "------------------------------------------------------");
                Console.WriteLine(indent + "------------------------------------------------------");*/
            }
            for(int i=0;i<nSteps;i++)
            //for (double x = a; x<=b-step; x+=step)
            {
                double x = a + step * i;
                x_plus_dx = x + step;
                if (x_plus_dx > b)
                    x_plus_dx = b;
                arguments[indexOfIntegratingVariable] = x;
                double func1 = func(arguments);
                arguments[indexOfIntegratingVariable] = x_plus_dx;
                //Integration xi
                double func2 = func(arguments);
                if (double.IsNaN(func1) || double.IsNaN(func2))
                    ;

                dS = step * (func2+func1) / 2.0;
                output += dS;
                if (double.IsNaN(output))
                    ;
                if (log)
                {
                   /* Console.WriteLine(indent + "--------------------Integration step----------------------");
                    Console.WriteLine(indent + "x = " + x);
                    Console.WriteLine(indent + "f(x) = " + func1);
                    Console.WriteLine(indent + "x+dx = " + x_plus_dx);
                    Console.WriteLine(indent + "f(x+dx) = " + func2);
                    Console.WriteLine(indent + "dS = " + dS);
                    Console.WriteLine(indent + "S =" + output);*/
                }
            }
            if (log)
            {
                //Console.WriteLine("RESULT = " + output);
            }
            if (double.IsNaN(output))
                output = output;
            return output;
        }
        static double XiFunction(double r,double component)
        {
            double integral = DefiniteIntegral(Xi_SubintegrativeExpression, -a, a, 0, new List<double>() { r, component });
            double output = integral / (2.0 * a);
            int y = 0;
           // if (double.IsNaN(output))
            //    output = 0;
            //return 0;
            return output;
        }
        static double Xi_SubintegrativeExpression(List<double> args)
        {
            double output = 0;
            double x = args[0];
            double r = args[1];
            int component = (int)args[2];
            for (int i = 0; i < numberOfComponents; i++)
            {
                if (double.IsNaN(Math.Sqrt(a * a - 2.0 * x * r + r * r)))
                    ;
                output += molarFractions[i] * LJpotential(Math.Sqrt(a*a-2.0*x*r+r*r), component, i);
                if (double.IsNaN(output))
                    ;
            }
            return z * output;
        }
        static double LJpotential(double r, int component1, int component2)
        {
            double sigma_r = sigmaMatrix[component1][component2] / r;

            double output= 4.0 * epsilonMatrix[component1][component2] * (Math.Pow(sigma_r,12.0) - Math.Pow(sigma_r, 6.0));
            if (double.IsNaN(output))
                ;
            if (output == 0)
                ;
            return output;
        }
        static StreamWriter xi_r_writer;
        static void GrapXiFromR()
        {
            return;
            string line = "";
            double step = 10E-12;
            line += a + ";";
            for(double r = 0; r < a; r += step)
            {
                double xi_r = XiFunction(r, 0);
                line += xi_r + ";";
            }
            xi_r_writer.WriteLine(line);
        }
        static double BisectionSolve(Func<double[],double> func,double accurasy, double a, double b, int xVariableIndexInArgs, List<double> additionalFuncArgs)
        {
            Console.WriteLine("Bisection soving start...");
            //Решение уравнений методом бисекции
                double maxIterations = 1000;
                double currentRoot = (a + b) / 2;
                int iter = 0;
                while (Math.Abs(a - b) > accurasy && maxIterations > 0)
                {
                    double[] aFuncArgs = new double[additionalFuncArgs.Count+1];
                    double[] bFuncArgs = new double[additionalFuncArgs.Count+1];
                    double[] middleFuncArgs = new double[additionalFuncArgs.Count+1];
                    int additionalArgsIndex = 0;
                    for (int i = 0; i < additionalFuncArgs.Count + 1; i++)
                    {
                        if (i == xVariableIndexInArgs)
                        {
                            aFuncArgs[i] = a;
                            bFuncArgs[i] = b;
                            middleFuncArgs[i] = (a + b) / 2;
                        }
                        else
                        {
                            aFuncArgs[i] = additionalFuncArgs[additionalArgsIndex];
                            bFuncArgs[i] = additionalFuncArgs[additionalArgsIndex];
                            middleFuncArgs[i] = additionalFuncArgs[additionalArgsIndex];
                            additionalArgsIndex++;
                        }
                    }
                    double f_a = func(aFuncArgs);
                    double f_b = func(bFuncArgs);
                    double f_middle = func(middleFuncArgs);
                    Console.WriteLine("Iteration: " + iter);
                    Console.WriteLine("f("+ aFuncArgs[1] +";" + aFuncArgs[2] + ";" + aFuncArgs[3] +") ="+ f_a);
                    Console.WriteLine("f(" + middleFuncArgs[1] + ";" + middleFuncArgs[2] + ";" + middleFuncArgs[3] + ") =" + f_middle);
                    Console.WriteLine("f(" + bFuncArgs[1] + ";" + bFuncArgs[2] + ";" + bFuncArgs[3] + ") =" + f_b);
                    maxIterations--;
                    iter++;
                if (double.IsNaN(f_a) && !double.IsNaN(f_b))
                {
                    /*a = (a + b) / 2;
                    currentRoot = (a + b) / 2;
                    continue;*/
                    f_a = 100000000000;
                }
                if (!double.IsNaN(f_a) && double.IsNaN(f_b))
                {
                    /*b= (a + b) / 2;
                    currentRoot = (a + b) / 2;
                    continue;*/
                    f_b = 100000000000;
                }


                    if (f_a * f_b > 0)
                    {
                    if (iter == 0)
                        return 0;
                    else
                        return (a + b) / 2;

                    }
                    if (f_a == 0)
                        return a;
                    if (f_b == 0)
                        return b;
                    if (f_a * f_middle > 0)
                        a = (a + b) / 2;
                    else
                        b = (a + b) / 2;
                    currentRoot = (a + b) / 2;
                }
                return currentRoot;
        }
    }
}