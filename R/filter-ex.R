df <- data.frame(A = c(0.5, 1.0))

df$Numlvl <- factor(df$A, levels = c(0.5, 1.0))

df$ManChar <- factor(df$A, levels = c("0.50", "1.0"))
