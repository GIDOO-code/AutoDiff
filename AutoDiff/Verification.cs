using System;
using System.Linq;
using static System.Console;
using static System.Math;

namespace AutoDiff{
    public  class VerificationMan {
        static private double eps=1.0e-10; 

        public VerificationMan( ){ }

        public void Verification_Functions( double p, double q, double[] ps, double[] qs ) {
            WriteLine( $"==================== Verification ==================== p:{p}  q:{q}");
            
            AD_Man ADM = new AD_Man();

            double PP = Math.PI/p, sg, pw;
            AD x1 = new AD(p.ToString());
            AD x2 = new AD(q.ToString());
            int N = ps.Length;
            double DN=N;
            AD[] x1s = new AD[N];
            AD[] x2s = new AD[N];

            AD y = +x1;
            _Verify( y, p, p, 1);
            
            y = -x1;
            _Verify( y, p, -p, -1);

            y = x1+x2;
            _Verify( y, p,q, p+q, 1, 1);
            
            y = x1-x2;
            _Verify( y, p,q, p-q, 1, -1);       
            y = x1*x2;
            _Verify( y, p,q, p*q, q, p);

            y = x1/x2;
            _Verify( y, p,q, p/q, 1/q, -p/(q*q));

            y = x1^q;
            _Verify( y, p,q, pw=Math.Pow(p,q), q*pw/p, pw*Math.Log(p) );
            
            WriteLine();
            y = AD.Exp(x1);
            _Verify( y, p, Math.Exp(p), Math.Exp(p) );

            //If the argument is negative, an ArgumentException will occur.
            // Partial differential value is invalid.
            y = AD.Sqrt(x1);
            _Verify( y, p, Math.Sqrt(p), 1/2.0/Math.Sqrt(p) );

            //If the argument is negative, an ArgumentException will occur.
            // Partial differential value is invalid.
            y = AD.Log(x1);
            _Verify( y, p, Math.Log(p), 1.0/p );

            //If the argument is negative, an ArgumentException will occur.
            // Partial differential value is invalid.
            y = AD.Log(x1,2.0);
            _Verify( y, p, Math.Log(p,2.0), 1.0/p/Math.Log(2.0) );

            y = AD.Pow(x1,x2);
            _Verify( y, p,q, pw=Math.Pow(p,q), q*pw/p, pw*Math.Log(p) );

            y = x1^x2;
            _Verify( y, p,q, pw=Math.Pow(p,q), q*pw/p, pw*Math.Log(p) );

            y = AD.Sin(x1);
            _Verify( y, PP, Math.Sin(PP), Math.Cos(PP) );

            y = AD.Cos(x1);
            _Verify( y, PP, Math.Cos(PP), -Math.Sin(PP) );

            y = AD.Tan(x1);
            _Verify( y, PP, Math.Tan(PP), 1.0/Math.Pow(Math.Cos(PP),2) );

            y = AD.Tanh(x1);
            _Verify( y, PP, Math.Tanh(PP), 1.0-Math.Pow(Math.Tanh(PP),2) );



            y = AD.Abs(x1);
            _Verify( y, p, Math.Abs(p), (p>0)? 1: -1 );

            y = AD.Max(x1,x2);
            _Verify( y, p,q, Math.Max(p,q), ((p>q)? 1: 0), ((p<q)? 1: 0) );

            y = AD.Min(x1,x2);
            _Verify( y, p,q, Math.Min(p,q), ((p>q)? 0: 1), ((p<q)? 0: 1) );

            y = AD.Sigmoid(p); 
            _Verify( y, p, sg=1/(1+Math.Exp(-p)), sg*(1-sg) );

            y = AD.ReLU(p); 
            _Verify( y, p, Math.Max(0,p), ((p<0)? 0: 1) );

            WriteLine();
            y = AD.Average(x1s); 
            var psL = ps.ToList();
            var qsL = qs.ToList();
            _Verify( y, ps, pw=psL.Sum()/N, Enumerable.Repeat<double>(1.0/DN,N).ToArray() );   

            y = AD.Sum(x1s); 
            _Verify( y, ps, psL.Sum(), Enumerable.Repeat<double>(1.0,N).ToArray() );       

            y = AD.InnerProd(x1s,x2s);
            double pd=0.0;
            for( int k=0; k<N; k++ ) pd += ps[k]*qs[k];
            _Verify( y, ps, qs, pd, qs, ps );     
        }



        private string TFCheck( double v, double vr ) => (Abs(v-vr)<=eps)? "ok": "▼";
        private void _Verify( AD y, double x1, double vr, double dr ) {
            try{ 
                AD_function ADF=y.ADF;
                double v = ADF.AD_func1(x1);
                double d = ((func_d)ADF.AD_dFuns[0])(x1);
                string st =  $" {y.opName.PadRight(6)} v:{v:f6} -- vr:{vr:f6} {TFCheck(v,vr)}";
                st += $"     d:{d:f6} -- dr:{dr:f6}     {TFCheck(d,dr)}" ;
                WriteLine(st);
            }
            catch( Exception e){ WriteLine(); WriteLine( "\r"+ $"{e.Message}\r{e.StackTrace}" ); }
        }

        private void _Verify( AD y, double x1,  double x2, double vr, double dr1, double dr2 ) {
            try{ 
            AD_function ADF=y.ADF;
            double v = ADF.AD_func2(x1,x2);
            double d1 = ((func_dd)ADF.AD_dFuns[0])(x1,x2);
            double d2 = ((func_dd)ADF.AD_dFuns[1])(x1,x2);
            string st =  $" {y.opName.PadRight(6)} v:{v:f6} -- vr:{vr:f6} {TFCheck(v,vr)}";
            st += $"     d:{d1:f6} -- dr:{dr1:f6}     {TFCheck(d1,dr1)}" ;
            st += $"     d:{d2:f6} -- dr:{dr2:f6}     {TFCheck(d2,dr2)}" ;
            WriteLine(st);
            }
            catch( Exception e){ WriteLine(); WriteLine( "\r"+ $"{e.Message}\r{e.StackTrace}" ); }
        }
        private void _Verify( AD y, double[] xs, double vr, double[] drs ) {
            try{ 
                AD_function ADF=y.ADF;
                double v = ADF.AD_funcN(xs);
                WriteLine( $" {y.opName.PadRight(6)} v:{v:f6} -- vr:{vr:f6} {TFCheck(v,vr)}" );
                for( int k=0; k<xs.Length; k++ ){
                    double d = ((func_dn)ADF.AD_dFuns[k])(xs);
                    string st = $"     d:{d:f6} -- dr:{drs[k]:f6}     {TFCheck(d,drs[k])}" ;
                    WriteLine(st);
                }  
            }
            catch( Exception e){ WriteLine(); WriteLine( "\r"+ $"{e.Message}\r{e.StackTrace}" ); }
        }
        private void _Verify( AD y, double[] x1s,  double[] x2s, double vr, double[] dr1s, double[] dr2s ){
            try{ 
                AD_function ADF=y.ADF;
                int N = x1s.Length;
                double v = ADF.AD_funcNN(x1s,x2s);
                WriteLine( $" {y.opName.PadRight(6)} v:{v:f6} -- vr:{vr:f6} {TFCheck(v,vr)}" );
                for( int k=0; k<N; k++ ){
                    double d = ((func_dndnk)ADF.AD_dFuns[k])(x1s,x2s,k);
                    string st = $"     d:{d:f6} -- dr:{dr1s[k]:f6}     {TFCheck(d,dr1s[k])}" ;
                    WriteLine(st);
                }
                for( int k=0; k<N; k++ ){
                    double d = ((func_dndnk)ADF.AD_dFuns[k+N])(x1s,x2s,k);
                    string st = $"     d:{d:f6} -- dr:{dr2s[k]:f6}     {TFCheck(d,dr2s[k])}" ;
                    WriteLine(st);
                }
            }
            catch( Exception e){ WriteLine(); WriteLine( "\r"+ $"{e.Message}\r{e.StackTrace}" ); }
        }

    }
}
 