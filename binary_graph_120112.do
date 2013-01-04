*clear all* Right now this is doing interactions.  Need to add squared later.* This is similar to the 'grinter' command from 
* http://myweb.uiowa.edu/fboehmke/methods.html
* but Boehmke does not provide proper marginal effects for probit models,
* He simply provides dy/dx assuming a linear relationship


*Define the program name
capture program drop binary_graphprogram define binary_graphversion 11args graph ci
set more off*Size of sample from multivariate normal for each loop iterationlocal j = 5000*Save the command line entry to use the names later in the graphscalar names = e(cmdline)local cmd = e(cmd)/*Find the number of RHS vars*/matrix beta = e(b)scalar columns = colsof(beta)*Calculate the median of all the RHS variables*Not sure exactly how this is working, but the appropriate*number of scalars will be produced belowloc indvars=e(indvars)                if("`indvars'"==".") {        loc indvars : colnames e(b)        loc indvars : list uniq indvars        loc constant _cons        loc indvars : list indvars - constant        }*Note this will produce a 1X(k+1) vector of medians called 'median'*The first element will be blank*The second element will be the variable we are taking the partial of*The third element will be the one the graph is modifying*The fourth element will be the median of the interaction term (meaningless)*The fifth through k+1 elements will be the medians of the control variables,*in the order they were placed in the regression*Also, the code is going to add in vectors for the min and max of the RHS varsmatrix RHS_median = .  matrix min_vec = .matrix max_vec = .                 foreach i of varlist `indvars' {        qui sum `i',d        scalar `i'_m=r(p50)        matrix RHS_median = (RHS_median, `i'_m)                 scalar `i'_l=r(min)        matrix min_vec = (min_vec, `i'_l)                scalar `i'_h=r(max)        matrix max_vec = (max_vec, `i'_h)             di `i'         }             *This extra step to convert the elements of RHS_median back to scalars is needed*because the above scalars are defined in terms of varlist names, not numbers*The second scalar will be the variable we are taking the partial of*The third scalar will be the one the graph is modifying*The fourth scalar will be the median of the interaction term (meaningless)*The fifth through k scalar will be the medians of the control variables,*in the order they were placed in the regressionlocal p = 2while `p' <= (columns) { scalar RHS_m`p' = RHS_median[1,`p']local p = `p'+1}* Scale will be used to make the graph
if ("`graph'" == "interaction"){       scalar min = min_vec[1,3]scalar max = max_vec[1,3]
}

if ("`graph'" == "square"){       
scalar min = min_vec[1,2]scalar max = max_vec[1,2]
}
scalar scale = max-min*     ****************************************************************  **       Take `j' draws from the estimated coefficient vector and     **       variance-covariance matrix.                                     **     ****************************************************************  *preservelocal b = columnsdrawnorm b_1-b_`b', n(`j') means(e(b)) cov(e(V)) clear*     ****************************************************************  **       To calculate the desired quantities of interest we need to set  **       up a loop.  This is what we do here.                            **     ****************************************************************  **       First, specify what you quantities should be saved and what     **       these quantities should be called.                              **     ****************************************************************  *postutil clearpostfile mypost  diff_hat diff_lo diff_hi using sim , replace            *     ****************************************************************  **       Start loop.  Let `a' be the modifying variable Z and let this   **       run from min to max in the desired increments.                  **     ****************************************************************  *                                  /*Set the size of the matrix for stata */set mat `j'/*Create a matrix of the simulated values from the betas and call it 'X'*/mkmat b_1-b_`b', matrix(X)/*Begin the loop*/local a=min while `a' <= max + scale/500 { /*Generate a vector of the median of the x's*/if ("`graph'" == "interaction"){       *The first 3 elements in 'medians' are the median of x1, x2 = a, and the third is x1*a*RHS_m2 is the median of the partial variablematrix medians = (RHS_m2,`a',`a'*RHS_m2) *Now add the medians of all the controls to the vector 'medians'local e = 5while `e' <= (columns) { matrix medians = (medians,RHS_m`e') local e = `e'+1}}if ("`graph'" == "square"){       *The first 2 elements in 'medians' are the variable of interest and its square matrix medians = (`a',`a'*`a') *Now add the medians of all the controls to the vector 'medians'local e = 4while `e' <= (columns) { matrix medians = (medians,RHS_m`e') local e = `e'+1}
}









*Finally add a constant to the vector of mediansmatrix medians = (medians,1)/*Calculate the marginal effect*/*Multiply the betas by the median values of the X's	matrix x_b0 =  X*medians'*Save the above vector as variables.  x_betahat is a 1X`j' *vector of simulated betas times the X's	svmat x_b0, names(x_betahat)***********************************
* For Probit Model*Send to the normal density since the ME is * normalden(X'B)B for the probit model.  See Greene p. 775 (6th ed.)  if("`cmd'"=="probit") gen prob=normalden(x_betahat)    	***********************************
* For Logit Model
* The ME for the logit model is lamba(1-lambda)B. See Greene p. 775 (6th ed.)	

if("`cmd'"=="logit") {
	*gen exp_xb = exp(x_betahat)
	*gen lambda = exp_xb/(1+exp_xb)
	gen prob = (exp(x_betahat)/(1+exp(x_betahat)))*(1-(exp(x_betahat)/(1+exp(x_betahat))))
}






************************************ stuff is the marginal effect for a particular observation  


* For models with interaction terms
if ("`graph'" == "interaction") gen stuff = prob*(b_1 + b_3*`a')if ("`graph'" == "square") gen stuff = prob*(b_1 + b_2*2*`a')* Save the mean and 95 % CI for each value of a        egen diff=mean(stuff)        tempname  diff_hat diff_lo diff_hi   if ("`ci'" == "ci_90"){           _pctile stuff, p(5,95)		/*for a 90 percent ci*/}else{    _pctile stuff, p(2.5,97.5)		/*for a 95 percent ci*/}
    scalar `diff_lo'= r(r1)    scalar `diff_hi'= r(r2)      scalar `diff_hat'=diff * Post the values that have been generated           post mypost (`diff_hat') (`diff_lo') (`diff_hi')                   drop x_betahat  prob  diff  stuff 
    
    
    
    * We are going to do the loop 500 times to get a pretty graph    local a=`a'+ (scale/500)     display "." _c    } postclose mypost*     ****************************************************************  *                                  *       Call on posted quantities of interest to make the graph     **     ****************************************************************  *                                  use sim, cleargen MV = ((_n-1)/500)*scale + min
quietly drop if _n > 501

*Generate the names of the variables to be used in the graphlocal name3 =word(names,3)local name4 =word(names,4)* Generate the label for the graphlabel variable diff_hat "Marginal Effect"

if ("`ci'" == "ci_90"){       label variable diff_lo "90% Confidence Interval"}else{
label variable diff_lo "95% Confidence Interval"}if ("`graph'" == "interaction"){       graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) clcolor(black) ///        ||   line diff_lo  MV, clpattern(dash) clwidth(thin) clcolor(red) ///        ||   line diff_hi  MV, clpattern(dash) clwidth(thin) clcolor(red) ///        ||  ,   ///            yline(0) ///            legend(order(1 2)) ///            title("Marginal Effect of `name3' as `name4' Changes ") ///            xtitle("`name4'") ///            ytitle("Marginal Effect of `name3'") ///
            xlabel(minmax)
            
}

if ("`graph'" == "square"){       graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) clcolor(black) ///        ||   line diff_lo  MV, clpattern(dash) clwidth(thin) clcolor(red) ///        ||   line diff_hi  MV, clpattern(dash) clwidth(thin) clcolor(red) ///        ||  ,   ///            yline(0) ///            legend(order(1 2)) ///            title("Marginal Effect of `name3' ") ///            xtitle("`name3'") ///            ytitle("Marginal Effect ") ///
			xlabel(minmax)
}			
			            end









/*
* An example of how to use the program


*Load datawebuse union, clearprobit union age grade not_smsa 

*Generate the interaction termgen age_grade = age*grade
gen age_2 = age^2
gen grade_2 = grade^2*Run the probit*probit union age grade age_grade not_smsa 
probit union age age_2 grade not_smsa 

* For some reason I can't get grinter to work after probit* grinter weight , inter(weight_mpg) const02(mpg)

*binary_graph square ci_95
*binary_graph interaction ci_95binary_graph square ci_90
*binary_graph interaction ci_90


*/
