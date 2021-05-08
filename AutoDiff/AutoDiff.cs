using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using static System.Console;
using static System.Math;

// Top-down automatic differentiation.
// Automatic Differentiation of Parallel OpenMP Programs.
// Compile the expression and calculate the value and partial derivative.


namespace AutoDiff{

    #region AD_Man
    public class AD_Man{
        //AD_Man manages the structure of an expression from a bird's-eye view.
        //FuncLst and FuncRevLst are arranged in the order of the input / output relationship of the function.
        //Nodelst is an array of variables that determine the order of functions.

        private List<AD>            NodeLst = new List<AD>();
        private List<AD_function>   FuncLst;
        private List<AD_function>   FuncRevLst;
        private int _computeStep=0;

        public AD_Man( ){
            NodeLst = new List<AD>();
            AD.pNodeLst = NodeLst;
            AD._ID0_ = 0;
            _computeStep=0;
        }

      #region Compile,Compute,Differetial
        public void AD_Compile( bool checkB=false ){ // <<..... Optional .....
            if(_computeStep>=2)  return;
                    if(checkB)  NodeLst.ForEach(P=> WriteLine( $"nodeLst1 {P.ToStringSP()}") );

            // Semiorder between variables based on reference relationship.
            NodeLst.ForEach(p=>p._level=1);
            FuncLst = new List<AD_function>();
            NodeLst.ForEach(N=> {if(N.ADF!=null) FuncLst.Add(N.ADF); } );
            bool SetF=true;
            int  loopC=0;
            int maxLoop = NodeLst.Count * 100;
            while( SetF && ++loopC<maxLoop ){
                SetF=false;
                foreach( var F in FuncLst ){
                    int lvl= F.ADz._level;
                    F.ADPs.ForEach( q => { if(q._level<=lvl){q._level=lvl+1; SetF=true;} } );
                }
#if false
                foreach( var F in FuncLst ){
                    string st = $"loop:{loopC}  F.ADz.opName:{F.ADz.opName}  level:{F.ADz._level}"+ "\r\n";
                    F.ADPs.ForEach( q => st+= $"        {q.opName}  level:{q._level}"+ "\r\n" );
                    WriteLine(st);
                }
                WriteLine("-----------------------------------------------");
#endif
            }
            if( loopC>=maxLoop ) throw new ArithmeticException( "Expression is looping." );
            
            int LMin = FuncLst.Min( f=> f.ADz._level );
            FuncLst.ForEach( f=> f.ADz._level-=LMin );
            FuncLst.Sort( (f,g)=> -(f.ADz._level-g.ADz._level));
                
            FuncRevLst = FuncLst.ConvertAll(p=>p);
            FuncRevLst.Reverse();
                    if(checkB){
                        WriteLine();
                        NodeLst.ForEach(P=> WriteLine( $"nodeLst2 {P.ToStringSP()}") );
                    }
            _computeStep = 1;
            return;
        }

        public void AD_Compute( AD[] Vars, double[] Points, bool printSw=false ){
            if(_computeStep==0)  AD_Compile();
            
            { //--- forward ---
                if(printSw){ WriteLine(); WriteLine("\r*** Forward mode ****"); }
                // Values ​​can be obtained by calculating in the order of expression definitions.
                for( int k=0; k<Vars.Length; k++ )  Vars[k].Val = Points[k];    //Set value for variables
                NodeLst.ForEach( P => P.Dif=0.0 );
                FuncLst.ForEach(F => F.ComputeFunc(printSw:printSw) );          //Function value calculation
            }

            { //--- Reverse ---
                if(printSw){ WriteLine(); WriteLine("\r*** Reverse mode ****"); }
                var F0 = FuncRevLst[0];
                NodeLst.ForEach( P => P.Dif=0.0 );
                F0.ADz.Dif= 1.0;
                FuncRevLst.ForEach(F => F.ComputeDifFunc(printSw) );        //Calculation of partial differential value
            }

            _computeStep = 2;
            return;
        }


      #endregion Compile,Compute,Differetial

        public void check_nodeLst( string ttl=""){
            WriteLine();
            WriteLine( $"\r\r*** NodeLst ***  ttl:{ttl}");
            NodeLst.ForEach( p=> WriteLine( $" \t {p.ToStringSP(p._level)}" ) );
        }
    }
    #endregion AD_Man

    #region delegate
    public delegate double func_d(  double x1 );                            //Func(double,double)
    public delegate double func_dd( double x1, double x2 );                 //Func(double,double,double)
    public delegate double func_dn( double[] xs );                          //Func(double[],double)
    public delegate double func_dndn( double[] xs, double[] ys );           //Func(double[],double[],double)
    public delegate double func_dndnk( double[] xs, double[] ys, int k );   //Func(double[],double[],int,double)
    #endregion delegate

    #region AD_function
    public class AD_function{
        static public int   _ID0_=0;
        private int         ID=-1;
        public  VerificationMan VeriMan=null;
        public  List<AD>    ADPs;                 // Input variables
        public  AD          ADz;                  // Output variable

        public func_d      AD_func1  = null;       //1 variable function
        public func_dd     AD_func2  = null;       //2-variable function
        public func_dn     AD_funcN  = null;       //N-variable function
        public func_dndn   AD_funcNN = null;       //2-Vector function

        public List<object> AD_dFuns = new List<object>();        //Partial differential functions


        #region Constructor
        public AD_function(){
            ID=_ID0_++;
            ADPs = new List<AD>();
        }

        // ----- Definition(function, output variable, input variable, partial differential function) -----
        public AD_function( func_d AD_func1, AD z, AD w, func_d df1 ): this(){
            this.AD_func1 = AD_func1;
            ADz = z;
            ADPs.Add(w);
            AD_dFuns.Add(df1);
        }
        public AD_function( func_d AD_func1, AD z, AD w, func_dd df2 ): this(){
            this.AD_func1 = AD_func1;
            ADz = z;
            ADPs.Add(w);
            AD_dFuns.Add(df2);
        }
        public AD_function( func_dd AD_func2, AD z, AD u, AD w, func_dd df2a, func_dd df2b ): this(){
            this.AD_func2 = AD_func2;
            ADz = z;
            ADPs.Add(u); ADPs.Add(w);       //not change order
            AD_dFuns.Add(df2a); AD_dFuns.Add(df2b);

        }
        public AD_function( func_dn ADfuncN, AD z, AD[] ws, List<func_dn> dfN ): this(){
            this.AD_funcN = ADfuncN;
            ADz = z;
            ADPs.AddRange(ws);
            AD_dFuns.AddRange(dfN);
        }

        public AD_function( func_dndn ADfuncNN, AD z, AD[] us, AD[] ws, List<func_dndnk> dfdfN ): this(){
            this.AD_funcNN = ADfuncNN;
            ADz = z;
            ADPs.AddRange(us);
            ADPs.AddRange(ws);
            AD_dFuns.AddRange(dfdfN); 
        }
        #endregion Constructor

        public string ToStringSP( int lvl=0 ){ 
            string stSP = ADz.ToStringSP(ADz._level);
            AD_dFuns.ForEach( df => WriteLine( df ) );
            return stSP;
        }




        // ----- Function, calculation, execution -----
        public void ComputeFunc( bool printSw=false ){
            if(AD_func1!=null){
                AD ADP=ADPs[0];
                ADz.Val = AD_func1(ADP.Val);
                if(printSw) WriteLine( $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = func1( ADp.Val(#{ADP.ID})[{ADP.Val:f6}]" );
            }
            else if(AD_func2!=null){
                AD ADP0=ADPs[0];
                AD ADP1=ADPs[1];
                ADz.Val = AD_func2( ADP0.Val, ADP1.Val );
                if(printSw){
                    string st = $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = ";
                    st += $"func2( ADp0.Val(#{ADP0.ID})[{ADP0.Val:f6}], ADp1.Val(#{ADP1.ID})[{ADP1.Val:f6}] )";
                    WriteLine(st);
                }
            }

            else if(AD_funcN!=null){
                double[] Q = ADPs.ConvertAll(q=>q.Val).ToArray();
                ADz.Val = AD_funcN( Q );
                if(printSw){
                    string st = "";
                    ADPs.ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                    WriteLine( $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = func2( ADPs.Val:{st} )" );
                }
            }
            else if(AD_funcNN!=null){
                var Q = ADPs.ConvertAll(q=>q.Val);
                int n = Q.Count/2;
                var Q0 = 
                ADz.Val = AD_funcNN( Q.Take(n).ToArray(), Q.Skip(n).ToArray() );
                if(printSw){
                    string st = "";
                    ADPs.Take(n).ToList().ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                    st += " / ";
                    ADPs.Skip(n).ToList().ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                    WriteLine( $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = func2( ADPs.Val:{st} )" );
                }
            }
        }

        public void ComputeDifFunc( bool printSw=false ){
            if(AD_func1!=null){                         //1 variable function                    
                AD ADP=ADPs[0];
                var FD0 = AD_dFuns[0];
                if( FD0 is func_d ){
                    func_d FD = FD0 as func_d;
                    ADP.Dif = ADz.Dif*FD(ADP.Val);
                    if(printSw){
                        string st =  $"\tF#{ID} ADp.Dif(#{ADP.ID})[{ADP.Dif:f6}] = ";
                        st += $"ADz.Dif(#{ADz.ID})[{ADz.Dif:f6}] X Dfunc1( ADp.Val(#{ADP.ID})[{ADP.Val:f6}] )";
                        WriteLine(st);
                    }
                }
                else{
                    throw new InvalidOperationException( $" Invalid definition of expression {ADz.opName}." );
                }
            }
            else if(AD_func2!=null){                    //2 variable function
                AD ADP0=ADPs[0];
                AD ADP1=ADPs[1];
                for( int k=0; k<2; k++ ){
                    var FDk = AD_dFuns[k];
                    if( FDk is func_dd ){
                        func_dd FD = FDk as func_dd;
                        ADPs[k].Dif += ADz.Dif*FD(ADP0.Val,ADP1.Val);
                        if(printSw){
                            string st = $"\tF#{ID} ADps{k}.Dif(#{ADPs[k].ID})[{ADPs[k].Dif:f6}] =";
                            st += $" ADz.Dif(#{ADz.ID})[{ADz.Dif:f6}] X ";
                            st += $"Dfunc2( ADPs0.Val(#{ADP0.ID})[{ADP0.Val:f6}], ADPs1.Val(#{ADP1.ID})[{ADP1.Val:f6}] )";
                            WriteLine(st); 
                        }
                    }
                    else{
                        throw new InvalidOperationException( $" Invalid definition of expression {ADz.opName}." );
                    }
                }
            }
            else if(AD_funcN!=null){
                double[] Q = ADPs.ConvertAll(q=>q.Val).ToArray();
                ADz.Val = AD_funcN( Q );
                if(printSw){
                    int n = Q.Length;
                    for( int k=0; k<n; k++ ){               //N variable function
                        string st = "";
                        ADPs.ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                        WriteLine( $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = funcn( ADPs.Val:{st} )" );
                    }
                }
            }
            else if(AD_funcNN!=null){                       //2-Vector function
                var Q = ADPs.ConvertAll(q=>q.Val);
                int n = Q.Count/2;
                ADz.Val = AD_funcNN( Q.Take(n).ToArray(), Q.Skip(n).ToArray() );
                if(printSw){
                    string st = "";
                    ADPs.Take(n).ToList().ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                    st += " / ";
                    ADPs.Skip(n).ToList().ForEach(P => { st += $", (#{P.ID})[{P.Val:f6}]"; } );
                    WriteLine( $"\tV#{ID} ADz.Val(#{ADz.ID})[{ADz.Val:f6}] = funcNN( ADPs.Val:{st} )" );
                }
            }
            else{
                throw new InvalidOperationException( $" Invalid definition of expression {ADz.opName}." );
            }
       
        }

    }
    #endregion AD_function
   
    #region AD
    public class AD{
        private const double _delta = 1e-13;
        static public int    _ID0_=0;
        static public List<AD> pNodeLst;

        public  int          ID=-1;
        public　int          _level;
        public  string       opName="var";
        public  AD_function  ADF;
        public  double       Val;                // Function value
        public  double       Dif;                // Partial differential value
        public  AD_Man       pADM=null;
        private AD[]         pVars=null;

        #region Constructor
        public AD(){
            ID=_ID0_++;
            opName = $"v{ID}";
            if(pNodeLst==null) throw new NullReferenceException( $"NullReferenceException: {opName}.ADM" );
            pNodeLst.Add(this);
            pADM=null;
        }

        public AD( double w ): this(){
            opName  =  w.ToString();
            Val     = w;
        }

        public void Evaluate( AD[] Vars, double[] Points, bool printSw=false ){
            this.pVars = Vars;
            if(ADF==null) throw new NullReferenceException( $"NullReferenceException: {opName}.ADM" );
            else   pADM.AD_Compute( Vars, Points, printSw );
        }


        public AD( string expression ): this(){
            opName = expression;
        }
        #endregion Constructor 



      #region Auxiliary routine        
        public double GetDifferential( ){
            return  Dif;
        }              

        public List<(double,double)> GetValDiffLst( AD[] Vars=null ){
            AD[] qVars = Vars?? pVars;
            List<(double,double)> ValDiffLst = (qVars.ToList()).ConvertAll(p=>(p.Val,p.Dif)).ToList();
            return ValDiffLst;
        }
        public double[] GetValLst( AD[] Vars=null ){
            AD[] qVars = Vars?? pVars;
            double[] ValLst = (qVars.ToList()).ConvertAll(p=>p.Val).ToArray();
            return ValLst;
        }
        public double[] GetDiffLst( AD[] Vars=null ){
            AD[] qVars = Vars?? pVars;
            double[] diffLst = (qVars.ToList()).ConvertAll(p=>p.Dif).ToArray();
            return diffLst;
        }

        public string ToStringSP( int lvl=0 ){ 
            string stSP = new String(' ',_level);
            stSP += $"ID:{ID}-L{_level}  opName:{opName.ToString().PadLeft(8)}  Val:{Val:f6}  Dif:{Dif:f6}";
            return stSP;
        }
      #endregion Auxiliary routine

      #region Unary operator, Conversion operator   
        static public implicit operator AD( double v ){// 暗黙の型変換
            var z = new AD(v);
            return z;
        }     

        static public AD operator +( AD x ){     // + Conversion operator
            var z = new AD("op+");
            z.ADF = new AD_function( (w)=>w, z, x, (w)=>1 );
            return z;
        }
        
        static public AD operator -( AD x ){    // - Conversion operator         
            var z = new AD("op-");
            z.ADF = new AD_function( (w)=>-w, z, x, (w)=>-1 );
            return z;
        }
      #endregion Unary operator, Conversion operator   
         
      #region Binary operator
        static public AD operator +( AD x, AD y ){      // + operator overload
            var z = new AD("op+op");
            z.ADF = new AD_function( (u,w)=>u+w, z, x, y, (u,w)=>1, (u,w)=>1 );
            return z;
        }　         
        static public AD operator -( AD x, AD y ){      // - operator overload
            var z = new AD("op-op");
            z.ADF = new AD_function( (u,w)=>u-w, z, x, y, (u,w)=>1, (u,w)=>-1 );
            return z;
        }
        static public AD operator *( AD x, AD y ){      // * operator overload
            var z = new AD("op*op");
            z.ADF = new AD_function( (u,w)=>u*w, z, x, y, (u,w)=>w, (u,w)=>u );
            return z;
        }       
        static public AD operator /( AD x, AD y ){      // / operator overload
            var z = new AD("op/op");
            z.ADF = new AD_function( (u,w)=>u/w, z, x, y, (u,w)=>1/w, (u,w)=>-u/(w*w) );
            return z;
        }

        // In AD, ^ is used for exponentiation. Enclose the target term in ().
        static public AD operator ^( AD x, AD y ){      // ^ operator overload    //------------???
            return Pow(x,y);
        }
        static public AD operator ^( AD x, double p ){      // ^ operator overload    //------------???
            AD x2 = new AD(p);
            return Pow(x,x2);
        }

        #endregion Binary operator

      #region Math functions       
        static public AD Exp( AD x ){       // Exp
            var z = new AD("Exp");
            z.ADF = new AD_function( (w)=>Math.Exp(w), z, x, (w)=>Math.Exp(w) );
            return z;
        }        
        static public AD Sqrt( AD x ){      // Sqrt
            var z = new AD("Sqrt");
#if aDEBUG
            func_d _fSqrt = (w) => {
                if(w<=0) throw new ArgumentException( $" Invalid Argument {z.opName}." );
                return Math.Sqrt(w);
            };
            z.ADF = new AD_function( _fSqrt, z, x, (w)=>0.5/Math.Sqrt(w) );
#else
            z.ADF = new AD_function( (w)=>Math.Sqrt(w), z, x, (w)=>0.5/Math.Sqrt(w) );
#endif
            return z;
        }      
        static public AD Log( AD x ){       // Log           
            var z = new AD("Log");
#if aDEBUG
            func_d _fLog = (w) => {
                if(w<=0) throw new ArgumentException( $" Invalid Argument {z.opName}." );
                return Math.Log(w);
            };
            z.ADF = new AD_function( _fLog, z, x, (w)=>1.0/w );
#else
            z.ADF = new AD_function( (w)=>Math.Log(w), z, x, (w)=>1.0/w );
#endif
            return z;
        }      
        static public AD Log( AD x, double a ){// Log(x,a)
            var z = new AD("Log_a");
#if aDEBUG
            func_d _fLoga = (w) => {
                if(w<=0||a<=0) throw new ArgumentException( $" Invalid Argument {z.opName}." );
                return Math.Log(w)/Math.Log(a);
            };
            z.ADF = new AD_function(_fLoga, z, x, (w)=>1/w/Math.Log(a));
#else
            z.ADF = new AD_function( (w)=>Math.Log(w)/Math.Log(a), z, x, (w)=>1/w/Math.Log(a));
#endif
            return z;
        }
        static public AD Pow( AD x, AD y ){ //  x^y         
            var z = new AD("Pow");
#if aDEBUG
            func_dd _fPow_dw = (u,w) => {
                if(w<=0) throw new ArgumentException( $" Invalid Argument {z.opName}." );
                return Math.Pow(u,w)*Math.Log(u);
            };
            z.ADF = new AD_function( (u,w)=>Math.Pow(u,w), z, x, y, (u,w)=>w*Math.Pow(u,w-1), _fPow_dw );
#else
            z.ADF = new AD_function( (u,w)=>Math.Pow(u,w), z, x, y, (u,w)=>w*Math.Pow(u,w-1), (u,w)=>Math.Pow(u,w)*Math.Log(u) );
#endif
            return z;
        }
        static public AD Sin( AD x ){       // Sin
            var z = new AD("Sin");
            z.ADF = new AD_function( (w)=>Math.Sin(w), z, x, (w)=>Math.Cos(w) );
            return z;
        }
        static public AD Cos( AD x ){       // Cos
            var z = new AD("Cos");
            z.ADF = new AD_function( (w)=>Math.Cos(w), z, x, (w)=>-Math.Sin(w) );
            return z;
        }
        static public AD Tan( AD x ){       // Tan
            var z = new AD("Tan");           
            z.ADF = new AD_function( (w)=>Math.Tan(w), z, x, (w)=> {
                double cos = Math.Cos(w);
                return 1/(cos*cos); } );
            return z;
        }        
        static public AD Tanh( AD x ){      // Tanh
            var z = new AD("Tanh");
            func_d  _tanh = (w)=> { double wt=Math.Tanh(w); return 1.0- wt*wt; };
            z.ADF = new AD_function( (w)=>Math.Tanh(w), z, x, _tanh );
            return z;
        }      
        static public AD Abs( AD x ){       // Abs
            var z = new AD("Abs");
            z.ADF = new AD_function( (w)=>Math.Abs(w), z, x, (w)=>(w<0? -1: 1) );
            return z;
        }     
        static public AD Max( AD x, AD y ){ // Max
            AD z = new AD("Max");
            z.ADF = new AD_function( (u,w)=>Math.Max(u,w), z, x, y, (u,w)=>(u>w?1:0), (u,w)=> (u>w?0:1) );
            return z;
        }        
        static public AD Min( AD x, AD y ){ // Min関数
            AD z = new AD("Min");
            z.ADF = new AD_function( (u,w)=>Math.Min(u,w), z, x, y, (u,w)=>(u<w?1:0), (u,w)=>(u<w?0:1) );
            return z;
        }     
        static public AD Sigmoid( AD x ){   // Sigmoid    
            var z = new AD("Sigmoid");
            func_d  _sig = (w)=> { double zz=1/(1+Math.Exp(-w)); return zz*(1.0-zz); };
            z.ADF = new AD_function( (w)=>1/(1+Math.Exp(-w)), z, x, _sig );
            return z;
        }
        static public AD ReLU( AD x ){      // ReLU            
            var z = new AD("ReLU");
            z.ADF = new AD_function( (w)=>Math.Max(0,w), z, x, (w)=>((w<=0)? 0:1) );
            return z;
        }
        static public AD Sum( AD[] xs ){    // Sum           
            var z = new AD("Sum");
            int N=xs.Length;
            List<func_dn> _dfN = Enumerable.Repeat<func_dn>((ws)=>1,N).ToList();
            z.ADF = new AD_function( (ws)=>ws.ToList().Sum(), z, xs, _dfN );
            return z;
        }    
        static public AD Average( AD[] xs ){ // Average          
            var z = new AD("Average");
            int N=xs.Length;
            double _dn = 1.0/N;
            List<func_dn> dfN = Enumerable.Repeat<func_dn>((ws)=>_dn,N).ToList();
            z.ADF = new AD_function( (ws)=>ws.ToList().Average(), z, xs, dfN );
            return z;
        }

        static public AD InnerProd( AD[] xs, AD[] ys ){// InnerProd
            var z = new AD("InnerProd");
            if(xs.Length!=ys.Length) throw new ArgumentException( $"ArgumentException: {z.opName}.ADM" );
            var N = Math.Min(xs.Length, ys.Length);

            func_dndn  _fnn = (us,ws)=> {
                double tt=0;
                for( int k=0; k<N; k++ ) tt += us[k]*ws[k];
                return tt;
            };
            List<func_dndnk> dfnnL = new List<func_dndnk>();

            for( int k=0; k<N; k++ )  dfnnL.Add( (us,ws,k)=>ws[k] );
            for( int k=0; k<N; k++ )  dfnnL.Add( (us,ws,k)=>us[k] );

            z.ADF = new AD_function( _fnn, z, xs, ys, dfnnL );
            return z;
        }
        #endregion Math functions

    }
    #endregion AD
}
