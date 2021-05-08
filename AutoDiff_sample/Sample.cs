using System;
using System.Linq;
using static System.Console;

using AutoDiff;

namespace AutoDiff {
    
    // *==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*
    //  Sample
    // *==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*


    public  class SampleMan {

        public SampleMan( ){ }
        public  void Sample1(){ 
            WriteLine( "==================== AD_Sample1 ====================" );
            AD_Man ADM = new AD_Man();          // AutoDiff manager. Defined before starting AD

            AD x1=new AD( );
            AD x2=new AD( );

            // Attention!   (x^n) : adjust the priority of operators.
            AD y = x1 * AD.Exp( -0.5 * ( (x1^2) + (x2^2) ) );                       // <-- #####-1

            //Set after function definition. Necessary! 
            y.pADM = ADM;




            AD[]   Vars   = { x1, x2 };
            double[] Points = {  2, 0.5};
            
            //Compile and calculate functions and partial derivatives.
            y.Evaluate( Vars, Points );                                             // <-- #####-2
                   //AD_Man.check_nodeLst( "AD_Man.AD_Compute  --- for debug" );
          
            

                // [type-1] How to refer to values ​​and partial derivatives.
                    string st = $" type-1  y.Val:{y.Val:f6}";
                    st += $"  x1.Val:{x1.Val:f6}  x2.Val:{x2.Val:f6}";
                    st += $"  x1.Dif:{x1.Dif:f6}  x2.Dif:{x2.Dif:f6}";
                    WriteLine(st);

                // [type-2] How to refer to values ​​and partial derivatives.
                var VDlst = y.GetValDiffLst( Vars:Vars );     
                    string st2 = $" type-2  y.Val:{y.Val:f6}";                  
                    foreach( var (P,k) in VDlst.Select((v,i)=>(v,i)))  st2+=$"  x{k+1}.Val:{P.Item1:f6}";
                    foreach( var (P,k) in VDlst.Select((v,i)=>(v,i)))  st2+=$"  x{k+1}.Dif:{P.Item2:f6}";
                    WriteLine(st2);
            //-------------------------------------------




            // The expression is defined by #####-1.
            // Calculate function and partial derivatives using compiled expressions(#####-2).
            for( double d=1; d<=5; d++){ 
                Points[0] = d;
                y.Evaluate( Vars, Points );
                
                var VDlst2 = y.GetValDiffLst( );        //If the argument is null, use the "Evaluate" argument.
                    string st3 = $" d:{d}  y.Val:{y.Val:f6}";                  
                    foreach( var (P,k) in VDlst2.Select((v,i)=>(v,i)))  st3+=$"  x{k+1}.Val:{P.Item1:f6}";
                    foreach( var (P,k) in VDlst2.Select((v,i)=>(v,i)))  st3+=$"  x{k+1}.Dif:{P.Item2:f6}";
                    WriteLine( st3 );
            }
        }   





      #region NewtonRaphson
        public void NewtonRaphson_1(){
            WriteLine( "==================== NewtonRaphson_1 ... y = 5.0-(x^2) ====================" );
            AD_Man ADM = new AD_Man();          // AD manager (Declare first)
 
            AD x = new AD();            
            // In AD, ^ is used for exponentiation. Enclose the target term in ().
            AD y = 5.0-(x^2);                   //(eq)  y=5.0-x*x or y=5.0-AD.Pow(x,2);
            y.pADM = ADM;                       // <- necessary!

            double x0 = _NewtonRaphson( y, x, 1.0, 20, 1.0e-20 );
            WriteLine( $"   solution:{x0:f10}  ... Sqrt(5.0)  ");
        }


        public void NewtonRaphson_2(){
            WriteLine( "==================== NewtonRaphson_2 ... y = (x^2) - 4*x + 3 ====================" );
            AD_Man ADM = new AD_Man();          // AD manager (Declare first)

            AD x = new AD();                
            AD y = (x^2) - 4*x + 3;             //(eq) y = AD.Pow(x,2)  - 4*x + 3;         
            y.pADM = ADM;                       // <- necessary!

            double x1 = _NewtonRaphson( y, x, 0, 20, 1.0e-20 );
            double x2 = _NewtonRaphson( y, x, 7, 20, 1.0e-20 );
            WriteLine( $"   solution _1:{x1:f10}  _2:{x2:f10}");
        }


        private double _NewtonRaphson( AD y, AD x, double initGuess, int maxIterations=20, double eps=1.0e-10 ){
            AD[]   Vars   = {x};
            double[] Points = new double[1];
                      
            double Xpre, Xsol=initGuess;
            int k=0;
            do{ 
                Xpre = Xsol;
                Points[0] = Xsol;
                y.Evaluate( Vars, Points, printSw:false ); //Compile and calculate functions and partial derivatives.
                Xsol -= y.Val / x.Dif;

                WriteLine( $"   loop:{k} Xsol:{Xsol:f10} ");
            }while( Math.Abs(Xpre-Xsol)>eps &&  ++k<maxIterations );

            return Xsol;
        }
      #endregion NewtonRaphson


    }

}
 