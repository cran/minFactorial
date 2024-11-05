
minFactDesign<-function(Levels_of_the_Factors,lower_bound,upper_bound){

  vectors <- lapply(Levels_of_the_Factors, gen.level)

  generate_recursive <- function(vectors, current_combination = NULL, index = 1) {

    if (index > length(vectors)) {

      return(matrix(unlist(current_combination), ncol = length(vectors), byrow = TRUE))

    } else {

      current_vector <- vectors[[index]]

      new_combinations <- lapply(current_vector, function(x) {

        new_current_combination <- append(current_combination, list(x))

        generate_recursive(vectors, new_current_combination, index + 1)
      })
      return(do.call(rbind, new_combinations))
    }
  }


  initial_design <- generate_recursive(vectors)

  generate_permutations <- function(numbers, subset_size) {

    perms <- vector("list")

    permute <- function(current_permutation, remaining_numbers) {

      if (length(current_permutation) == subset_size) {

        perms <<- append(perms, list(current_permutation))

      } else {

        for (i in seq_along(remaining_numbers)) {

          next_permutation <- c(current_permutation, remaining_numbers[i])

          next_remaining_numbers <- remaining_numbers[-i]

          permute(next_permutation, next_remaining_numbers)
        }
      }
    }

    permute(c(), numbers)

    return(perms)
  }

  permutation<-generate_permutations(1:nrow(initial_design),nrow(initial_design))

  permutation_of_rows <- do.call(rbind, permutation)

  Total_Minimally_Changed_Factorial_Run_Orders <- 0

  Minimally_Changed_Factorial_Run_Orders<-list()

  All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor<- list()

  D_values <- c()

  Dt_values <- c()

  Max_D_value <- -Inf

  Max_Dt_value <- -Inf

  D_optimal_designs <- list()

  Dt_optimal_designs <- list()

  Max_Trend_factor_value <- -Inf

  Number_of_Designs_Max_Trend_Factor <- 0

  Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range <- list()

  for (j in 1:nrow(permutation_of_rows)) {

    design <- initial_design[permutation_of_rows[j, ], , drop = FALSE]

    count <- numeric(ncol(design))

    for (k in 1:ncol(design)) {

      for (l in 2:nrow(design)) {

        if (design[l - 1, k] != design[l, k]) {

          count[k] <- count[k] + 1
        }
      }
    }

    Orthogonal_polynomial<-poly(1:nrow(initial_design),degree = 1,simple = TRUE)

    if (sum(count) == nrow(initial_design) - 1) {

      Total_Minimally_Changed_Factorial_Run_Orders<-Total_Minimally_Changed_Factorial_Run_Orders+1

      Minimally_Changed_Factorial_Run_Orders<-append(Minimally_Changed_Factorial_Run_Orders,list(design))


      mean_column <- matrix(1, nrow = nrow(initial_design), ncol = 1)

      design_matrix <- cbind(mean_column, design)

      A1<-t(design_matrix)%*%design_matrix

      D_value<- det(A1)

      D_values <- c(D_values, D_value)

      A2<-t(design_matrix)%*% Orthogonal_polynomial

      A3<-t(A2)

      A4<-t(Orthogonal_polynomial)%*%Orthogonal_polynomial

      Dt_value<- round(det(rbind(cbind(A1, A2),cbind(A3, A4))),2)

      threshold <- 1e-4

      if (Dt_value < threshold) {

        Dt_value <- 0
      }

      Dt_values <- c(Dt_values, Dt_value)

      Trend_factor <- round((Dt_value / D_value)^(1 / ncol(design_matrix)),2)


      All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor<-append(All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor,
                                                                                list(list(design,D_value=D_value,Dt_value=Dt_value,Trend_factor=Trend_factor)))

      if (D_value > Max_D_value) {


        Max_D_value <- D_value

        D_optimal_designs <- list(list(design, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor))

      } else if (D_value == Max_D_value) {


        D_optimal_designs <- append(D_optimal_designs, list(list(design, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor)))
      }

      if (Dt_value > Max_Dt_value) {

        Max_Dt_value <- Dt_value

        Dt_optimal_designs <- list(list(design, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor))

      } else if (Dt_value == Max_Dt_value) {

        Dt_optimal_designs <- append(Dt_optimal_designs, list(list(design, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor)))
      }

      if (!is.nan(Trend_factor) && Trend_factor > Max_Trend_factor_value) {

        Max_Trend_factor_value <- Trend_factor

        Number_of_Designs_Max_Trend_Factor <- 1

      } else if (!is.nan(Trend_factor) && Trend_factor == Max_Trend_factor_value) {

        Number_of_Designs_Max_Trend_Factor <- Number_of_Designs_Max_Trend_Factor + 1
      }

      if (!is.nan(Trend_factor) && Trend_factor >= lower_bound && Trend_factor <= upper_bound) { # trend factor should not be NAN


        design_with_Trend_factor <- list(design,D_value=D_value,

                                         Dt_value=Dt_value, Trend_factor = Trend_factor)

        Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range <- append(Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range, list(design_with_Trend_factor))
      }

    }

  }
  if (length(Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range) == 0) {

    Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range <- "No designs found within the specified range of trend factor"
  }


  if (length(unique(D_values)) == 1) {
    D_optimal_designs <- "No D optimal designs found."
  }


  if (length(unique(Dt_values)) == 1) {
    Dt_optimal_designs <- "No Dt optimal designs found."
  }

  return(list(Max_D_value = Max_D_value,
              D_optimal_designs=D_optimal_designs,
              Max_Dt_value = Max_Dt_value,
              Dt_optimal_designs = Dt_optimal_designs,
              Max_Trend_factor_value=Max_Trend_factor_value,
              Number_of_Designs_Max_Trend_Factor=Number_of_Designs_Max_Trend_Factor,
              Minimally_Changed_Factorial_Run_Orders=Minimally_Changed_Factorial_Run_Orders,
              Total_Minimally_Changed_Factorial_Run_Orders=Total_Minimally_Changed_Factorial_Run_Orders,
              All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor=All_Minimally_Changed_Factorial_Run_Orders_with_D_Dt_Trend_Factor,
              Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range = Minimally_Changed_Factorial_Run_Orders_in_trend_factor_range))

}

