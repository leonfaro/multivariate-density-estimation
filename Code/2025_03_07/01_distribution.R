# 01_distribution.R

distLibrary <- list(
  
  Normal = list(
    paramNames = c("mean","sd"),
    rFun = function(n, mean, sd) {
      if(length(mean)==1) mean <- rep(mean,n)  # replicate if scalar
      if(length(sd)==1)   sd   <- rep(sd,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- rnorm(1, mean=mean[i], sd=sd[i])  # single draws
      }
      out
    },
    dFun = function(x, mean, sd) {
      n <- length(x)
      if(length(mean)==1) mean <- rep(mean,n)
      if(length(sd)==1)   sd   <- rep(sd,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- dnorm(x[i], mean=mean[i], sd=sd[i])
      }
      out
    }
  ),
  
  Exponential = list(
    paramNames = c("rate"),
    rFun = function(n, rate) {
      if(length(rate)==1) rate <- rep(rate,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- rexp(1, rate=rate[i])
      }
      out
    },
    dFun = function(x, rate) {
      n <- length(x)
      if(length(rate)==1) rate <- rep(rate,n)
      out <- numeric(n)
      for(i in seq_along(x)) {
        out[i] <- dexp(x[i], rate=rate[i])
      }
      out
    }
  ),
  
  Maxwell = list(
    paramNames = c("sigma"),
    rFun = function(n, sigma) {
      if(length(sigma)==1) sigma <- rep(sigma,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- sigma[i]*sqrt(rchisq(1, df=3))
      }
      out
    },
    dFun = function(x, sigma) {
      n <- length(x)
      if(length(sigma)==1) sigma <- rep(sigma,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        if(x[i]>0) {
          out[i] <- sqrt(2/pi)*(x[i]^2/(sigma[i]^3))*exp(-x[i]^2/(2*sigma[i]^2))
        } else {
          out[i] <- 0
        }
      }
      out
    }
  ),
  
  Logistic = list(
    paramNames = c("mu","s"),
    rFun = function(n, mu, s) {
      if(length(mu)==1) mu <- rep(mu,n)
      if(length(s)==1)  s  <- rep(s,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- rlogis(1, location=mu[i], scale=s[i])
      }
      out
    },
    dFun = function(x, mu, s) {
      n <- length(x)
      if(length(mu)==1) mu <- rep(mu,n)
      if(length(s)==1)  s  <- rep(s,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- dlogis(x[i], location=mu[i], scale=s[i])
      }
      out
    }
  ),
  
  Lognormal = list(
    paramNames = c("meanlog","sdlog"),
    rFun = function(n, meanlog, sdlog) {
      if(length(meanlog)==1) meanlog <- rep(meanlog,n)
      if(length(sdlog)==1)   sdlog   <- rep(sdlog,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- rlnorm(1, meanlog=meanlog[i], sdlog=sdlog[i])
      }
      out
    },
    dFun = function(x, meanlog, sdlog) {
      n <- length(x)
      if(length(meanlog)==1) meanlog <- rep(meanlog,n)
      if(length(sdlog)==1)   sdlog   <- rep(sdlog,n)
      out <- numeric(n)
      for(i in seq_len(n)) {
        out[i] <- dlnorm(x[i], meanlog=meanlog[i], sdlog=sdlog[i])
      }
      out
    }
  )
)

create_config_entry <- function(distName, paramFun) {
  list(distName=distName, paramFun=paramFun)
}

config <- vector("list",25)

config[[1]] <- create_config_entry(
  "Normal",
  function() {
    list(mean=0, sd=1)
  }
)

config[[2]] <- create_config_entry(
  "Exponential",
  function(Y1) {
    r_ <- 1 + 0.05*Y1
    r_[r_<0.01] <- 0.01
    list(rate=r_)
  }
)

config[[3]] <- create_config_entry(
  "Maxwell",
  function(Y1,Y2) {
    s_ <- 1 + 0.01*Y2
    s_[s_<0.01] <- 0.01
    list(sigma=s_)
  }
)

config[[4]] <- create_config_entry(
  "Logistic",
  function(Y1,Y2,Y3) {
    mu_ <- 0.2*Y1 + 0.3*Y2 + 0.1*Y3
    s_  <- 1 + 0.05*abs(Y1+Y2)
    s_[s_<0.01] <- 0.01
    list(mu=mu_, s=s_)
  }
)

config[[5]] <- create_config_entry(
  "Normal",
  function(Y1,Y2,Y3,Y4) {
    m_  <- Y1 + 0.5*Y3
    sd_ <- 1 + 0.01*(Y2+Y4)
    sd_[sd_<0.01] <- 0.01
    list(mean=m_, sd=sd_)
  }
)

config[[6]] <- create_config_entry(
  "Lognormal",
  function(Y1,Y2,Y3,Y4,Y5) {
    ml_ <- 0.1*(Y1+Y2) + 0.01*Y5
    sd_ <- 0.5 + 0.001*abs(Y3+Y4)
    sd_[sd_<0.01] <- 0.01
    list(meanlog=ml_, sdlog=sd_)
  }
)

config[[7]] <- create_config_entry(
  "Exponential",
  function(Y1,Y2,Y3,Y4,Y5,Y6) {
    r_ <- 1 + 0.02*(Y1+Y4)
    r_[r_<0.01] <- 0.01
    list(rate=r_)
  }
)

config[[8]] <- create_config_entry(
  "Maxwell",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7) {
    s_ <- 1 + 0.05*Y7
    s_[s_<0.01] <- 0.01
    list(sigma=s_)
  }
)

config[[9]] <- create_config_entry(
  "Normal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8) {
    m_  <- 0.1*(Y2+Y3+Y6)
    sd_ <- 1 + 0.01*abs(Y5)
    sd_[sd_<0.01] <- 0.01
    list(mean=m_, sd=sd_)
  }
)

config[[10]] <- create_config_entry(
  "Logistic",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9) {
    mu_ <- 0.5*Y8
    s_  <- 1 + 0.01*abs(Y1+Y9)
    s_[s_<0.01] <- 0.01
    list(mu=mu_, s=s_)
  }
)

config[[11]] <- create_config_entry(
  "Lognormal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10) {
    ml_ <- 0.3 + 0.01*Y10
    sd_ <- 0.4 + 0.02*abs(Y2+Y3)
    sd_[sd_<0.01] <- 0.01
    list(meanlog=ml_, sdlog=sd_)
  }
)

config[[12]] <- create_config_entry(
  "Maxwell",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11) {
    s_ <- 1 + 0.01*Y6
    s_[s_<0.01] <- 0.01
    list(sigma=s_)
  }
)

config[[13]] <- create_config_entry(
  "Exponential",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12) {
    r_ <- 0.5 + 0.001*(Y9+Y11)
    r_[r_<0.01] <- 0.01
    list(rate=r_)
  }
)

config[[14]] <- create_config_entry(
  "Normal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13) {
    m_  <- Y12 + 0.2*Y13
    sd_ <- 1 + 0.01*abs(Y9+Y10)
    sd_[sd_<0.01] <- 0.01
    list(mean=m_, sd=sd_)
  }
)

config[[15]] <- create_config_entry(
  "Lognormal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14) {
    ml_ <- 0.1*Y14 + 0.001*Y8
    sd_ <- 0.3 + 0.002*abs(Y10)
    sd_[sd_<0.01] <- 0.01
    list(meanlog=ml_, sdlog=sd_)
  }
)

config[[16]] <- create_config_entry(
  "Logistic",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y15) {
    mu_ <- 0.2*Y5 + 0.3*Y15
    s_  <- 1 + 0.01*abs(Y14)
    s_[s_<0.01] <- 0.01
    list(mu=mu_, s=s_)
  }
)

config[[17]] <- create_config_entry(
  "Maxwell",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y15,Y16) {
    s_ <- 1 + 0.01*Y16
    s_[s_<0.01] <- 0.01
    list(sigma=s_)
  }
)

config[[18]] <- create_config_entry(
  "Normal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17) {
    m_  <- 0.5*Y17
    sd_ <- 1 + 0.02*abs(Y13)
    sd_[sd_<0.01] <- 0.01
    list(mean=m_, sd=sd_)
  }
)

config[[19]] <- create_config_entry(
  "Exponential",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18) {
    r_ <- 0.2 + 0.01*Y18
    r_[r_<0.01] <- 0.01
    list(rate=r_)
  }
)

config[[20]] <- create_config_entry(
  "Logistic",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19) {
    mu_ <- 0.3*Y19
    s_  <- 1 + 0.01*abs(Y8+Y18)
    s_[s_<0.01] <- 0.01
    list(mu=mu_, s=s_)
  }
)

config[[21]] <- create_config_entry(
  "Lognormal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20) {
    ml_ <- 0.2 + 0.01*(Y16+Y20)
    sd_ <- 0.3 + 0.02*abs(Y4)
    sd_[sd_<0.01] <- 0.01
    list(meanlog=ml_, sdlog=sd_)
  }
)

config[[22]] <- create_config_entry(
  "Maxwell",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20,Y21) {
    s_ <- 1 + 0.01*Y21
    s_[s_<0.01] <- 0.01
    list(sigma=s_)
  }
)

config[[23]] <- create_config_entry(
  "Normal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20,
           Y21,Y22) {
    m_  <- Y22 + 0.01*Y15
    sd_ <- 1 + 0.01*abs(Y12+Y14)
    sd_[sd_<0.01] <- 0.01
    list(mean=m_, sd=sd_)
  }
)

config[[24]] <- create_config_entry(
  "Exponential",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20,
           Y21,Y22,Y23) {
    r_ <- 0.1 + 0.01*Y23
    r_[r_<0.01] <- 0.01
    list(rate=r_)
  }
)

config[[25]] <- create_config_entry(
  "Lognormal",
  function(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,
           Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18,Y19,Y20,
           Y21,Y22,Y23,Y24) {
    ml_ <- 0.1*(Y23+Y24)
    sd_ <- 0.2 + 0.01*abs(Y18)
    sd_[sd_<0.01] <- 0.01
    list(meanlog=ml_, sdlog=sd_)
  }
)




