## docs/make_header.R -- draws docs/header.png (README banner)
suppressPackageStartupMessages(library(ggplot2))
set.seed(7)
K <- 6; p <- 14
alpha <- c(0.9, -0.7, 0.5, 0, 0.6, 0, -0.4, 0, 0, 0.3, 0, 0, -0.5, 0)
dev <- matrix(0, K, p)
for (j in which(alpha != 0)) { ks <- sample(K, sample(0:2, 1)); dev[ks, j] <- rnorm(length(ks), 0, 0.45) }
d <- expand.grid(k = 1:K, j = 1:p); d$theta <- alpha[d$j] + dev[cbind(d$k, d$j)]
d$type <- ifelse(dev[cbind(d$k, d$j)] != 0, "dev", ifelse(alpha[d$j] != 0, "shared", "zero"))
d$x <- d$j + (d$k - 3.5) * 0.09
g <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7A2E2E") +
  annotate("segment", x = 0.3, xend = 14.6, y = 0, yend = 0, colour = "white", alpha = 0.35) +
  geom_segment(data = data.frame(j = 1:p, a = alpha), aes(x = j - 0.35, xend = j + 0.35, y = a, yend = a),
               colour = "#F4ECEC", linewidth = 1.4) +
  geom_point(data = d, aes(x, theta, colour = type), size = 2.2) +
  scale_colour_manual(values = c(dev = "#F4B942", shared = "#F4ECEC", zero = "#C99C9C"), guide = "none") +
  annotate("text", x = 15.4, y = 0.55, label = "hiermetacox", hjust = 0, colour = "white",
           size = 13, fontface = "bold", family = "sans") +
  annotate("text", x = 15.45, y = 0.05, hjust = 0, colour = "#F4ECEC", size = 4.6, lineheight = 1.1,
           label = "Hierarchical nonconvex penalized stratified Cox models\nfor individual-participant-data meta-analysis") +
  annotate("text", x = 15.45, y = -0.55, hjust = 0, colour = "#F4B942", size = 3.8,
           label = "theta[kj] == alpha[j] + epsilon[kj]", parse = TRUE) +
  coord_cartesian(xlim = c(0.3, 30), ylim = c(-1.2, 1.3), expand = FALSE) + theme_void()
ggsave("docs/header.png", g, width = 12, height = 2.6, dpi = 150)
