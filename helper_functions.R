# Function to install a package if not installed and load
use_package <- function(p) {                                                                      
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
library(p, character.only = TRUE)
}

pif <- function (dgm, h = 0.2, n = 50, lims = c(0, 1, 0, 1)) 
{
    maxdim <- max(dgm$diag[, 1])
    pif_list <- list()
    for (i in 1:(maxdim + 1)) {
        list_label <- paste("dim_", (i - 1), sep = "")
        if (sum(dgm$diag[, 1] == (i - 1)) > 1) {
            tmp_pif <- MASS::kde2d((dgm$diag[(dgm$diag[, 1] == 
                (i - 1)), ])[, 2], (dgm$diag[(dgm$diag[, 1] == 
                (i - 1)), ])[, 3], h = h, n = n, lims = lims)
            tmp_pif$z[lower.tri(tmp_pif$z)] <- 0
            pif_list[[list_label]] <- tmp_pif$z/sum(tmp_pif$z)
        }
        else {
            tmp_pif <- MASS::kde2d((dgm$diag[(dgm$diag[, 1] == 
                (i - 1)), ])[2], (dgm$diag[(dgm$diag[, 1] == 
                (i - 1)), ])[3], h = h, n = n, lims = lims)
            tmp_pif$z[lower.tri(tmp_pif$z)] <- 0
            pif_list[[list_label]] <- tmp_pif$z/sum(tmp_pif$z)
        }
    }
    return(pif_list)
}

