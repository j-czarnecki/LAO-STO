      H_normal(1, 1) = -mu
      H_normal(2, 1) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(3, 1) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(4, 1) = -t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(-cmplx(0,1)*(0.86602540378443865d0*
     @ k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(5, 1) = 2*cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*
     @ sin(0.86602540378443865d0*k_x)
      H_normal(6, 1) = -t_rashba*(-exp(-cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(7, 1) = 0
      H_normal(8, 1) = -1.0d0/3.0d0*sqrt(6.0d0)*lambda_soc
      H_normal(9, 1) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),-1.0d0
     @ /2.0d0*sqrt(2.0d0)))
      H_normal(10, 1) = 0
      H_normal(11, 1) = 0
      H_normal(12, 1) = 0
      H_normal(1, 2) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(2, 2) = -mu
      H_normal(3, 2) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(4, 2) = -2*cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*
     @ sin(0.86602540378443865d0*k_x)
      H_normal(5, 2) = -t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(-cmplx(0,1)*(-0.86602540378443865d0*
     @ k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(6, 2) = -t_rashba*(-exp(-cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(7, 2) = (1.0d0/3.0d0)*sqrt(6.0d0)*lambda_soc
      H_normal(8, 2) = 0
      H_normal(9, 2) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(10, 2) = 0
      H_normal(11, 2) = 0
      H_normal(12, 2) = 0
      H_normal(1, 3) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(2, 3) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(3, 3) = -mu
      H_normal(4, 3) = t_rashba*(-exp(-cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(5, 3) = t_rashba*(-exp(-cmplx(0,1)*(0.86602540378443865d0
     @ *k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(6, 3) = -t_pi*exp(-cmplx(0,1)*k_y) - 2*t_sigma*exp(0.5d0*
     @ cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
      H_normal(7, 3) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(8, 3) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),(1.0d0
     @ /2.0d0)*sqrt(2.0d0)))
      H_normal(9, 3) = 0
      H_normal(10, 3) = 0
      H_normal(11, 3) = 0
      H_normal(12, 3) = 0
      H_normal(1, 4) = -t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(2, 4) = 2*cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*
     @ sin(0.86602540378443865d0*k_x)
      H_normal(3, 4) = t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(4, 4) = -mu
      H_normal(5, 4) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(6, 4) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(7, 4) = 0
      H_normal(8, 4) = 0
      H_normal(9, 4) = 0
      H_normal(10, 4) = 0
      H_normal(11, 4) = -1.0d0/3.0d0*sqrt(6.0d0)*lambda_soc
      H_normal(12, 4) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(1, 5) = -2*cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)
     @ *sin(0.86602540378443865d0*k_x)
      H_normal(2, 5) = -t_pi*exp(cmplx(0,1)*(0.86602540378443865d0*k_x -
     @ 0.5d0*k_y)) - t_sigma*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(3, 5) = t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(4, 5) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(5, 5) = -mu
      H_normal(6, 5) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(7, 5) = 0
      H_normal(8, 5) = 0
      H_normal(9, 5) = 0
      H_normal(10, 5) = (1.0d0/3.0d0)*sqrt(6.0d0)*lambda_soc
      H_normal(11, 5) = 0
      H_normal(12, 5) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(1, 6) = -t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(2, 6) = -t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(3, 6) = -t_pi*exp(cmplx(0,1)*k_y) - 2*t_sigma*exp(-0.5d0*
     @ cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
      H_normal(4, 6) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(5, 6) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(6, 6) = -mu
      H_normal(7, 6) = 0
      H_normal(8, 6) = 0
      H_normal(9, 6) = 0
      H_normal(10, 6) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(11, 6) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(12, 6) = 0
      H_normal(1, 7) = 0
      H_normal(2, 7) = (1.0d0/3.0d0)*sqrt(6.0d0)*lambda_soc
      H_normal(3, 7) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(4, 7) = 0
      H_normal(5, 7) = 0
      H_normal(6, 7) = 0
      H_normal(7, 7) = -mu
      H_normal(8, 7) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(9, 7) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(10, 7) = -t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0*
     @ k_x - 0.5d0*k_y)) - t_sigma*(exp(-cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(11, 7) = 2*cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*
     @ sin(0.86602540378443865d0*k_x)
      H_normal(12, 7) = -t_rashba*(-exp(-cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(1, 8) = -1.0d0/3.0d0*sqrt(6.0d0)*lambda_soc
      H_normal(2, 8) = 0
      H_normal(3, 8) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),-1.0d0
     @ /2.0d0*sqrt(2.0d0)))
      H_normal(4, 8) = 0
      H_normal(5, 8) = 0
      H_normal(6, 8) = 0
      H_normal(7, 8) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(8, 8) = -mu
      H_normal(9, 8) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(10, 8) = -2*cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)
     @ *sin(0.86602540378443865d0*k_x)
      H_normal(11, 8) = -t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(-cmplx(0,1)*(-0.86602540378443865d0*
     @ k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(12, 8) = -t_rashba*(-exp(-cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(1, 9) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),(1.0d0
     @ /2.0d0)*sqrt(2.0d0)))
      H_normal(2, 9) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(3, 9) = 0
      H_normal(4, 9) = 0
      H_normal(5, 9) = 0
      H_normal(6, 9) = 0
      H_normal(7, 9) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(3.0d0)
     @ *cmplx(0,1)*lambda_soc
      H_normal(8, 9) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(9, 9) = -mu
      H_normal(10, 9) = t_rashba*(-exp(-cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(11, 9) = t_rashba*(-exp(-cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
      H_normal(12, 9) = -t_pi*exp(-cmplx(0,1)*k_y) - 2*t_sigma*exp(0.5d0
     @ *cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
      H_normal(1, 10) = 0
      H_normal(2, 10) = 0
      H_normal(3, 10) = 0
      H_normal(4, 10) = 0
      H_normal(5, 10) = (1.0d0/3.0d0)*sqrt(6.0d0)*lambda_soc
      H_normal(6, 10) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(7, 10) = -t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(8, 10) = 2*cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)
     @ *sin(0.86602540378443865d0*k_x)
      H_normal(9, 10) = t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(10, 10) = -mu
      H_normal(11, 10) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(12, 10) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(1, 11) = 0
      H_normal(2, 11) = 0
      H_normal(3, 11) = 0
      H_normal(4, 11) = -1.0d0/3.0d0*sqrt(6.0d0)*lambda_soc
      H_normal(5, 11) = 0
      H_normal(6, 11) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),
     @ -1.0d0/2.0d0*sqrt(2.0d0)))
      H_normal(7, 11) = -2*cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y
     @ )*sin(0.86602540378443865d0*k_x)
      H_normal(8, 11) = -t_pi*exp(cmplx(0,1)*(0.86602540378443865d0*k_x
     @ - 0.5d0*k_y)) - t_sigma*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(
     @ -0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(9, 11) = t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(
     @ 0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(10, 11) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(11, 11) = -mu
      H_normal(12, 11) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(1, 12) = 0
      H_normal(2, 12) = 0
      H_normal(3, 12) = 0
      H_normal(4, 12) = lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(5, 12) = lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),(
     @ 1.0d0/2.0d0)*sqrt(2.0d0)))
      H_normal(6, 12) = 0
      H_normal(7, 12) = -t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*
     @ (-0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(8, 12) = -t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*
     @ (0.86602540378443865d0*k_x - 0.5d0*k_y)))
      H_normal(9, 12) = -t_pi*exp(cmplx(0,1)*k_y) - 2*t_sigma*exp(-0.5d0
     @ *cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
      H_normal(10, 12) = (1.0d0/2.0d0)*delta_tri - 1.0d0/3.0d0*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(11, 12) = (1.0d0/2.0d0)*delta_tri + (1.0d0/3.0d0)*sqrt(
     @ 3.0d0)*cmplx(0,1)*lambda_soc
      H_normal(12, 12) = -mu
