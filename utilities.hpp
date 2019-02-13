#include <string>
#include <vector>
#include <sstream>

#include "./muparser/include/muParser.h"

// Read an array of string split by a delim. 
// Return a new vector
std::vector<std::string>
split(const std::string & s, 
	  char delim
	  ) 
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) 
			elems.push_back(item);
    }
	return elems;
}

std::string expr="";
double wx=1, wy=1, Wx=1, Wy=1, H=1; 
// Function for source term in the inner domain
/*
double g_function(const bgeot::base_node & x){
    //Define p
    mu::Parser p;
    //Define variables
        // We define both coordinates x,y,z,
         //   * which are needed for the function,
         //   * variable you want. 
    double xx=x[0], yy=x[1], zz=x[2];
    p.DefineVar("x", &xx);
    p.DefineVar("y", &yy);
    p.DefineVar("z", &zz);
	p.DefineVar("wx", &wx);	
	p.DefineVar("wy", &wy);
	p.DefineVar("Wx", &Wx);
    p.DefineVar("Wy", &Wy);
    p.DefineVar("H", &H);
    //Define the expression
        // The string should be already modified
            // in the code 
    p.SetExpr(expr);
    return p.Eval();
}
*/
// Function for source term in the outer domain 
double f_function(const bgeot::base_node & x){

    //Define p
    mu::Parser p;
    //Define variables
        /* We define both coordinates x,y,z,
            * which are needed for the function.
            * Note that you can add in this file any other 
            * variable you want. */ 
    double xx=x[0], yy=x[1], zz=x[2];
    p.DefineVar("x", &xx);
    p.DefineVar("y", &yy);
    p.DefineVar("z", &zz);
	p.DefineVar("wx", &wx);	
	p.DefineVar("wy", &wy);
	p.DefineVar("Wx", &Wx);
    p.DefineVar("Wy", &Wy);
    p.DefineVar("H", &H);

    //Define the expression
        /* The string should be already modified
            * in the code */
    p.SetExpr(expr);
    return p.Eval();
}
