library("crs")
library("doParabar")
library("foreach")

source("boot_der.R")

df <- read.csv("contrast_face.csv")
vars <- colnames(df)
Y_names <- vars[10:369]
X <- df[c("BIS_Contrast", "Age", "Gender")]
# convert the two covariates to factors
X$Age = ordered(X$Age)
X$Gender = factor(X$Gender)

# Y = "R_V1_ROI"

for (Y in Y_names)
{
    f <- formula(paste(Y, "~ BIS_Contrast + ordered(Age) + factor(Gender)"))
    
    # automatic bandwidth selection, using least-squares CV.
    model_bw <- npglpreg(formula = f, data = df, bwtype = "auto", degree = 4, nmulti = 5, cv = "bandwidth")
    bws = model_bw$bws
    
    # first derivative and second derivative
    #model_first <- npglpreg(tydat = df[, Y], txdat = X, bwtype = "fixed", bws = bws, degree = 4, gradient.vec = 1)
    # model_second <- npglpreg(tydat = df[, Y], txdat = X, bwtype = "fixed", bws = bws, degree = 4, gradient.vec = 2)
    
    compute_mean_deriv(Y = df[, Y], X = X, bws = bws, degree = 4, gradient.vec = 1, B1 = 1000, B2 = 100)
    
    #first_deriv <- model_first$gradient
    #second_deriv <- model_second$gradient
    
    # hack: manual overwriting X and Y needed to workaround the bug of crs package
    # see plot.npglpreg in np.regression.glp.R for details.
    #model$x = X
    #model$y = df[, Y]
    
    plot(model, deriv = 0, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
    plot(model, deriv = 1, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
    plot(model, deriv = 2, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
}
