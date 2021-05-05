using System;
using System.Linq;
using static System.Console;

using AutoDiff;

namespace AutoDiff{
    public class Program {

        public delegate AD func_d( AD x1 );             //Func(double,double,double)
                  
        static void Main( string[] args ){

            // ==== Verification,  =====
            {
                VerificationMan Veri=new VerificationMan();

                Veri.Verification_Functions(5,2, new double[]{1,2,3,4}, new double[]{3,1,4,1} );
                Veri.Verification_Functions(-3,-2, new double[]{-1,2,-3,4}, new double[]{1,2,3,4} );
                Veri.Verification_Functions(3,10, new double[]{3,1,4,1}, new double[]{-1,2,-3,4} );
            }

            // ===== Sample =====
            { 
                SampleMan smpl = new SampleMan();
                smpl.Sample1();

                smpl.NewtonRaphson_1( );
                smpl.NewtonRaphson_2( );
            }
        }
    }
}
 