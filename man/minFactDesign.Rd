\name{minFactDesign}
\alias{minFactDesign}
\title{Construct all possible minimally changed factorial run orders for different experimental set ups}
\description{
The minFactDesign() function constructs all possible minimally changed factorial run orders for different experimental setups. The function takes inputs related to the levels of factors and the bounds for trend factor. Factor levels can be represented as integers, for example, -1 for low, 0 for medium, and 1 for high.User is expected to enter a vector of total number of levels for each factor to be considered in the experiment.If there are two factors, each with two levels, this should be entered as c(2,2). If there are three factors, each with two levels, this should be entered as c(2,2,2).If there are two factors, one with two levels and the other with three levels, it should be entered as c(2,3).The trend factor value lies between 0 to 1. Higher the values, lesser the influence of trend effects on the run order. The function then generates minimally changed factorial run orders and evaluates them based on D-optimality and Dt-optimality criteria. Additionally, it explores the impact of trend effects and identifies the designs that maximize the trend factor value.
}
\usage{
minFactDesign(Levels_of_the_Factors,lower_bound,upper_bound)
}
\arguments{
\item{Levels_of_the_Factors}{Number of Levels of Each Factor}
\item{lower_bound}{Lower bound of Trend Factor Value}
\item{upper_bound}{Upper bound of Trend Factor Value}
}
\value{
Returns a list with the following components:
\item{Total_Minimally_Changed_Factorial_Run_Orders}{Total number of all the minimally changed factorial run orders for the given inputs}
\item{Minimally_Changed_Factorial_Run_Orders}{All the minimally changed factorial run orders for the given inputs}
\item{All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor}{All the minimally changed factorial run orders with D, Dt and Trend Factor value for the given inputs}
\item{Max_D_value}{Maximum D-value within the generated minimally changed factorial run orders}
\item{D_optimal_designs}{Designs with the Maximum D-value within the generated minimally changed factorial run orders}
\item{Max_Dt_value}{Maximum Dt-value within the generated minimally changed factorial run orders}
\item{Dt_optimal_designs}{Designs with the Maximum Dt-value within the generated minimally changed factorial run orders}
\item{Max_Trend_factor_value}{Maximum Trend Factor Value for the generated minimally changed factorial run orders}
\item{Number_of_Designs_Max_Trend_Factor}{Number of minimally changed factorial run orders with maximum trend factor value}
\item{Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range}{Minimally changed factorial run orders within the specified range of trend factor}
}
\references{
Arpan Bhowmik, Eldho Varghese, Seema Jaggi and Cini Varghese (2015).Factorial experiments
with minimum changes in run sequences. Journal of the Indian Society of Agricultural Statistics,69(3), 243-255.

Arpan Bhowmik, Eldho Varghese, Seema Jaggi and Cini Varghese (2017).Minimally changed run sequences in factorial experiments. Communications in Statistics - Theory and Methods, 46(15), 7444-7459.

Arpan Bhowmik, Eldho Varghese, Seema Jaggi and Cini Varghese(2022).On the generation of factorial designs with minimum level changes. Communications in Statistics-Simulation and Computation, 51(6), 3400-3409.

Bijoy Chanda, Arpan Bhowmik, Seema Jaggi, Eldho Varghese, Anindita Datta, Cini Varghese, Namita Das Saha, Arti Bhatia and Bidisha Chakrabarti (2021). Minimal cost multifactor experiments for agricultural research involving hard-to-change factors. The Indian Journal of Agricultural Sciences, 91(7), 97-100.


Lieven Tack and Martina Vandebroek (2001). (Dt,C)-optimal run orders. Journal of Statistical Planning and Inference, 98, 293-310.


}
\examples{
# Two Factor each at Two Level
Result1 <- minFactDesign(c(2,2),0.8,0.9)

# Accessing results
Result1$Total_Minimally_Changed_Factorial_Run_Orders
Result1$Minimally_Changed_Factorial_Run_Orders
Result1$All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor
Result1$Max_D_value
Result1$D_optimal_designs
Result1$Max_Dt_value
Result1$Dt_optimal_designs
Result1$Max_Trend_factor_value
Result1$Number_of_Designs_Max_Trend_Factor
Result1$Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range

}

\keyword{Factorial Experiments}
