# Automatic differentiation Library 

(1) Reverse mode automatic differentiation.  
(2) Implemented compile function.  
    The calculation order of the function is determined from the input/output of the function.  
(3) Efficiently calculate function values ​​and partial differential values  

```
public  void Sample1(){ test0505
    AD_Man ADM = new AD_Man();          // AutoDiff manager. Defined before starting AD

    AD x1=new AD( );
    AD x2=new AD( );

    // Attention!   (x^n) : adjust the priority of operators.
    AD y = x1 * AD.Exp( -0.5 * ( (x1^2) + (x2^2) ) );

    //Set after function definition. Necessary! 
    y.pADM = ADM;

    AD[] Vars = { x1, x2 };
    double[] Points = { 2, 0.5};
    
    // Compile and calculate functions and partial derivatives.
    y.Evaluate( Vars, Points );
        
    // How to refer to values ​​and partial derivatives.
    string st = $" y.Val:{y.Val:f6}";
    st += $"  x1.Val:{x1.Val:f6}  x2.Val:{x2.Val:f6}";
    st += $"  x1.Dif:{x1.Dif:f6}  x2.Dif:{x2.Dif:f6}";
    WriteLine(st);
}
```
