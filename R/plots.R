#' Draws a plot for an MDP model.
#' @export
#' @title PlotMDPModel
#' @name PlotMDPModel
#' @param mdp.model MDP.model instance
#' @param separate Determines whether a separate plot for each
#' component is going to be produced
#' @return corresponding ggplot2 object
PlotMDPModel <- function (mdp.model, separate = FALSE) {
  y <- mdp.model$y
  s <- 0
  if (!is.null(mdp.model$theta$rho))
    s <- nrow(mdp.model$theta$rho)
  N <- ncol(y)
  T <- nrow(y) - s
  
  y.transposed <- data.frame(t(y))
  colnames(y.transposed) <- (-s+1):T
  y.df <- reshape::melt(cbind(id = as.factor(1:N), 
                              component = as.factor(mdp.model$components),
                              y.transposed), 
                        id = c('id', 'component'))
  colnames(y.df)[colnames(y.df)=='value'] <- 'y'
  colnames(y.df)[colnames(y.df)=='variable'] <- 't'
  
  
  if (!separate)
    return (ggplot2::ggplot(y.df, aes(x = t, y = y, 
                                      group = id, color = component)) +
      geom_point(size = 2) +
      geom_line(size = 1.2, alpha = 0.5))
  
  
  ggplot2::ggplot(y.df, aes(x = t, y = y, group = id, color = id)) +
    geom_point(size = 2) +
    geom_line(size = 1.2, alpha = 0.5) +
    facet_wrap(~component,
               nrow = 1) + guides(color=FALSE)
}