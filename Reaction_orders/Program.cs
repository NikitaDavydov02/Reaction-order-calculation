using System.Net.Http.Headers;
using System.Net.NetworkInformation;
using System.IO;
using System.Diagnostics.Tracing;

namespace Reaction_orders
{
    internal class Program
    {
       /* static double k = 1.38*Math.Pow(10, -23);//Boltzman constant
        static double K = Math.Pow(10, -10);
        static double T =15;//Temperature
        static double E0 = 400;//J/mol
        static double V = 0;//Volume
        static double h= 6.62 * Math.Pow(10, -34);
        static double Nav = 6.02 * Math.Pow(10, 23);
        static double R = 8.314;
        static double pi = Math.PI;
        static double z = 6;
        static double[] M;//Molar masses
        static double[] n;//Amounts
        static double[] gamma;
        static double[] chemPot;
        static double[] standartChemPot;
        static double[] Q;
        static double[][] LJparameters;
        static double[][] sigmaMatrix;
        static double[][] epsilonMatrix;
        static double[][] w;
        static double[][] beta;
        static double[][] noncorelatedPairsNumber;
        static double[][] corelatedPairsNumber;
        static int numberOfComponents = 3;
        static double[] m;
        static double[] N;
        static double[] molarFractions;
        static double[] chemPotForInitialPressure;
        static double a;
        static double Nsum;
        static double minV = 2;
        static double maxV = 3;
        static double nVSteps = 100;*/

        static bool log = false;
        static StreamWriter log_writer;
        static StreamWriter isotherm_writer;
        static StreamWriter chemicalPotential_writer;

        static double T = 15;
        static double P = 101325;
        
        delegate double Function(double x);
        
        static void Main(string[] args)
        {
            log_writer = new StreamWriter(File.Create("log.txt"));
            chemicalPotential_writer = new StreamWriter(File.Create("chemical_potential.txt"));

            //ReactionOrdersForDifferntCompositions();
            ReactionSystem sys = new ReactionSystem(2, 15, new double[3] { 12000, 12000, 0 });

            log_writer.Close();
            chemicalPotential_writer.Close();
            Console.WriteLine("Main is done");
         }
         static string isotherm_writer_line;
         /// <summary>
         /// Calculates reaction orders for different compositions of two component system
         /// </summary>
         /// <param name="steps"></param>
         /// <param name="n_sum"></param>
        static void ReactionOrdersForDifferntCompositions(int steps = 5, double n_sum=35000)
        {
            StreamWriter orderWriter = new StreamWriter("orders.csv");
            orderWriter.WriteLine("x0;x1;V;absolute;relative");
            double xStep = 0.98 / (double)steps;
            for(double x0 = 0.01; x0 < 1; x0 += xStep)
            {
                Console.WriteLine("-----------------------NEW COMPOSITION-------------------------");
               // V = FindVolumeCorrespondingToParticularPressure(101325);
               // UpdateParameters();
               // double p = CalculatePressureAtCurrentConditions();
               // UpdateEquilibriumComposition(n,chemPotForInitialPressure);
                //chemicalPotential_writer.WriteLine(x0+"     "+ gamma[0] + "    " + gamma[1] + "    " + gamma[2]+"    " + chemPot[0] + "    " + chemPot[1] + "    " + chemPot[2] + "    "+V);
               //  Console.WriteLine(x0 + "     " + gamma[0] + "    " + gamma[1] + "    " + gamma[2] + "    " + chemPot[0] + "    " + chemPot[1] + "    " + chemPot[2] + "    "+V);


                
                List<double> orders = CalculateReactionOrders(n_sum * x0, n_sum * (1 - x0));
                Console.WriteLine("Absolute order A: " + orders[0]);
                Console.WriteLine("Relative order A: " + orders[1]);
                Console.WriteLine("Absolute order B: " + orders[2]);
                Console.WriteLine("Relative order B: " + orders[3]);
                orderWriter.WriteLine(x0 + ";" + (1 - x0) + ";" + orders[0] + ";" + orders[1]+ ";" + orders[2] + ";" + orders[3] + ";");
                
            }
            orderWriter.Close();
        }
        /// <summary>
        /// Calculates all types of reaction orders for particular conditions
        /// </summary>
        /// <returns></returns>
         static List<double> CalculateReactionOrders(double n0, double n1)
         {
            List<double> output = new List<double>();
            double[] initial_composition = new double[3];
            initial_composition[0] = n0;
            initial_composition[1] = n1;
            initial_composition[2] = 0;

            ReactionSystem system = new ReactionSystem(P, T, initial_composition);
            double initialRate = system.CalculateCurrentReactionRate();
            double initialPressure = system.GetPressure();
            Console.WriteLine("Initial pressure: " + initialPressure);
            Console.WriteLine("Initial reaction rate:" + initialRate);
            double initV = system.V;
            

            //Absolute reaction order calculation
            initial_composition[0] = 1.01*n0;
            double dn_A = 0.01 * n0;
            system = new ReactionSystem(initV,T,initial_composition);
            double newPressure = system.GetPressure();
            double rateAfterNonIsobaricAddition = system.CalculateCurrentReactionRate();
            Console.WriteLine("Pressure after non-isobaric addition of A: " + newPressure);
            Console.WriteLine("Rate after non-isobaric addition of A: " + rateAfterNonIsobaricAddition);
            double V = system.V;
            double absoluteOrder = (Math.Log(rateAfterNonIsobaricAddition)-Math.Log(initialRate))/ (Math.Log(1.01 * n0 / V)- Math.Log(n0 / initV));
            Console.WriteLine("Absolute order of A: " + absoluteOrder);
            Console.WriteLine("Non-isobaric pressure for A: " + newPressure);

            //Relative order calculation
            system.SetPressure(P, FreedomDegree.N1);

            double pressureForRelativeOrderOfA = system.P;
            double rateAfterIsobaricAddition = system.CalculateCurrentReactionRate();
            V = system.V;
            Console.WriteLine("Reaction rate after isobaric addition: " + rateAfterIsobaricAddition);
            double relativeOrder = (Math.Log(rateAfterIsobaricAddition) - Math.Log(initialRate)) / (Math.Log(1.01 * n0 / V) - Math.Log(n0 / initV));
            Console.WriteLine("Relative order of A: " + relativeOrder);
            Console.WriteLine("Pressure for Relative order of A: " + pressureForRelativeOrderOfA);

            ///////////////////////////////////////
            //////////////////////////////////////
            ///////////////////////////////////////
            ///////////////////////////////////////
            //Absolute reaction order calculation of B
            initial_composition[0] = n0;
            initial_composition[1] = 1.01*n0;
            double dn_B = 0.01 * n1;
            system = new ReactionSystem(initV, T, initial_composition);
            double pressureNonIsobaricB = system.GetPressure();
            double rateAfterNonIsobaricAdditionofB = system.CalculateCurrentReactionRate();
            Console.WriteLine("Pressure after non-isobaric addition of B: " + pressureNonIsobaricB);
            Console.WriteLine("Rate after non-isobaric addition of B: " + rateAfterNonIsobaricAdditionofB);
            V = system.V;
            double absoluteOrderofB = (Math.Log(rateAfterNonIsobaricAdditionofB) - Math.Log(initialRate)) / (Math.Log(1.01 * n1 / V) - Math.Log(n1 / initV));
            Console.WriteLine("Absolute order of B: " + absoluteOrderofB);
            Console.WriteLine("Non-isobaric pressure for B: " + pressureNonIsobaricB);





            //Relative order calculation
            system.SetPressure(P, FreedomDegree.N0);
            double pressureForRelativeOrderOfB = system.P;
            double rateAfterIsobaricAdditionofB = system.CalculateCurrentReactionRate();
            V = system.V;
            Console.WriteLine("Reaction rate after isobaric addition of B: " + rateAfterIsobaricAdditionofB);
            double relativeOrderofB = (Math.Log(rateAfterIsobaricAdditionofB) - Math.Log(initialRate)) / (Math.Log(1.01 * n1 / V) - Math.Log(n1 / initV));
            Console.WriteLine("Relative order of B: " + rateAfterIsobaricAdditionofB);
            Console.WriteLine("Pressure for Relative order of B: " + pressureForRelativeOrderOfB);

            output.Add(absoluteOrder);
            output.Add(relativeOrder);
            output.Add(absoluteOrderofB);
            output.Add(relativeOrderofB);
            return output;
        }
         static void Isotherm(bool idealGas=false)
         {
            /* Console.WriteLine("Isotherm calculation start");
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
                     //if (idealGas)
                     //    F = -k * T * Calculate_LnZ_Gas();
                     //else
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
             Console.WriteLine("Isotherm calculation is done");*/
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
        
    }
}