using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Reaction_orders
{
    public static class MathFunctions
    {
        public static double BisectionSolve(Func<List<double>, double> func, double accurasy, double a, double b, int xVariableIndexInArgs, List<double> additionalFuncArgs)
        {
            Console.WriteLine("Bisection soving start...");
            //Решение уравнений методом бисекции
            double maxIterations = 1000;
            double currentRoot = (a + b) / 2;
            int iter = 0;
            while (Math.Abs(a - b) > accurasy && maxIterations > 0)
            {
                List<double> aFuncArgs = new List<double>();
                List<double> bFuncArgs = new List<double>();
                List<double> middleFuncArgs = new List<double>();
                int additionalArgsIndex = 0;
                for (int i = 0; i < additionalFuncArgs.Count + 1; i++)
                {
                    aFuncArgs.Add(0);
                    bFuncArgs.Add(0);
                    middleFuncArgs.Add(0);
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
                Console.WriteLine("f(" + aFuncArgs[1] + ";" + aFuncArgs[2] + ";" + aFuncArgs[3] + ") =" + f_a);
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
        public static double CalculateDerivative(Func<List<double>, double> func, double x, int xVariableIndexInArgs, List<double> additionalFuncArgs)
        {
            List<double> args = new List<double>();
            int additionalArgsIndex = 0;
            for (int i = 0; i < additionalFuncArgs.Count + 1; i++)
            {
                if (i == xVariableIndexInArgs)
                    args[i] = x;
                else
                {
                    args[i] = additionalFuncArgs[additionalArgsIndex];
                    additionalArgsIndex++;
                }
            }
            double dx = x * 1.01;
            args[xVariableIndexInArgs] = x - 2 * dx;
            double f1 = func(args);
            args[xVariableIndexInArgs] = x - dx;
            double f2 = func(args);
            args[xVariableIndexInArgs] = x + 2 * dx;
            double f3 = func(args);
            args[xVariableIndexInArgs] = x + dx;
            double f4 = func(args);
            double output = (8 * f4 - 8 * f2 - f3 + f1) / (12 * dx);
            return output;
        }
        public static double DefiniteIntegral(Func<List<double>, double> func, double a, double b, int indexOfIntegratingVariable, List<double> nonIntegratingVariablesValues, double nSteps = 1000)
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
            for (int i = 0; i <= nonIntegratingVariablesValues.Count; i++)
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
            for (int i = 0; i < nSteps; i++)
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

                dS = step * (func2 + func1) / 2.0;
                output += dS;
                if (double.IsNaN(output))
                    ;
                //if (log)
                //{
                    /* Console.WriteLine(indent + "--------------------Integration step----------------------");
                     Console.WriteLine(indent + "x = " + x);
                     Console.WriteLine(indent + "f(x) = " + func1);
                     Console.WriteLine(indent + "x+dx = " + x_plus_dx);
                     Console.WriteLine(indent + "f(x+dx) = " + func2);
                     Console.WriteLine(indent + "dS = " + dS);
                     Console.WriteLine(indent + "S =" + output);*/
                //}
            }
            //if (log)
            //{
                //Console.WriteLine("RESULT = " + output);
            //}
            if (double.IsNaN(output))
                output = output;
            return output;
        }
        public static double StirlingFormula(double N)
        {
            if (N == 0)
                return 1;
            return (N * Math.Log(N) - N);
        }
    }
}
