H_A1(1, 1) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(2, 1) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(3, 1) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(4, 1) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(5, 1) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(6, 1) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(7, 1) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_A1(8, 1) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(9, 1) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(10, 1) = 0
H_A1(11, 1) = 0
H_A1(12, 1) = 0
H_A1(13, 1) = 0
H_A1(14, 1) = 0
H_A1(15, 1) = 0
H_A1(16, 1) = (1.0d0/2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(cmplx(0,1)*( &
      -exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(9)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) - cmplx(0, &
      1)*exp(-cmplx(0,1)*k_y)*conjg(xi(10)) - cmplx(0,1)*exp(-cmplx(0,1 &
      )*k_y)*conjg(xi(13)) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(7 &
      ))
H_A1(17, 1) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (7)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0 &
      /3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(12)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx &
      (0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(10)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(11))
H_A1(18, 1) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      *conjg(xi(10)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg &
      (xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(6)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(13)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(19, 1) = 0
H_A1(20, 1) = 0
H_A1(21, 1) = 0
H_A1(22, 1) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(19)) + exp( &
      -cmplx(0,1)*k_y)*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(18)) - exp(-cmplx(0,1)*k_y)*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(23, 1) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(2)) + (1.0d0/2.0d0 &
      )*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)* &
      k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4))
H_A1(24, 1) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(19)) - exp( &
      -cmplx(0,1)*k_y)*conjg(xi(1)) + exp(-cmplx(0,1)*k_y)*conjg(xi(15 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(1, 2) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_A1(2, 2) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(3, 2) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(4, 2) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(5, 2) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(6, 2) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(7, 2) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(8, 2) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_A1(9, 2) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(10, 2) = 0
H_A1(11, 2) = 0
H_A1(12, 2) = 0
H_A1(13, 2) = 0
H_A1(14, 2) = 0
H_A1(15, 2) = 0
H_A1(16, 2) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (7)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0 &
      /3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(12)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx &
      (0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(10)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(11))
H_A1(17, 2) = (1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(9)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(11)) + (1.0d0/2.0d0 &
      )*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(14)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(8))
H_A1(18, 2) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg( &
      xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(7)) + (1.0d0/2.0d0) &
      *(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(10))
H_A1(19, 2) = 0
H_A1(20, 2) = 0
H_A1(21, 2) = 0
H_A1(22, 2) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(2)) + (1.0d0/2.0d0 &
      )*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)* &
      k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4))
H_A1(23, 2) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(1)) + (1.0d0/2.0d0 &
      )*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)* &
      k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(15)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(4))
H_A1(24, 2) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(1, 3) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(2, 3) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(3, 3) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(4, 3) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(5, 3) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(6, 3) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*k_y) - t_sigma*exp(0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(7, 3) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(8, 3) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(9, 3) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_A1(10, 3) = 0
H_A1(11, 3) = 0
H_A1(12, 3) = 0
H_A1(13, 3) = 0
H_A1(14, 3) = 0
H_A1(15, 3) = 0
H_A1(16, 3) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      *conjg(xi(10)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg &
      (xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(6)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(13)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(17, 3) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg( &
      xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(7)) + (1.0d0/2.0d0) &
      *(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(10))
H_A1(18, 3) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y))*conjg(xi(11)) + (1.0d0/2.0d0)*(-2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(14)) + (1.0d0/2.0d0)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(8)) + (1.0d0/2.0d0)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) &
      *conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(9))
H_A1(19, 3) = 0
H_A1(20, 3) = 0
H_A1(21, 3) = 0
H_A1(22, 3) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(19)) - exp( &
      -cmplx(0,1)*k_y)*conjg(xi(1)) + exp(-cmplx(0,1)*k_y)*conjg(xi(15 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(23, 3) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(24, 3) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*k_y))*conjg(xi(3)) + (1.0d0/2.0d0)*(-exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(15)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(4))
H_A1(1, 4) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(2, 4) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(3, 4) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 4) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(5, 4) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(6, 4) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(7, 4) = 0
H_A1(8, 4) = 0
H_A1(9, 4) = 0
H_A1(10, 4) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(11, 4) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(12, 4) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(13, 4) = (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-cmplx(0,1)*( &
      -exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(9)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + cmplx(0, &
      1)*exp(cmplx(0,1)*k_y)*conjg(xi(10)) + cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(13)) + cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(7))
H_A1(14, 4) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(11))
H_A1(15, 4) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      *conjg(xi(10)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg &
      (xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(6)) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(13)) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(16, 4) = 0
H_A1(17, 4) = 0
H_A1(18, 4) = 0
H_A1(19, 4) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      cmplx(0,1)*k_y)*conjg(xi(16)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(18)) - exp(cmplx(0,1)*k_y)*conjg(xi(2)) - 1.0d0/2.0d0* &
      exp(cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(20, 4) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(21, 4) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      cmplx(0,1)*k_y)*conjg(xi(1)) - exp(cmplx(0,1)*k_y)*conjg(xi(15)) &
      - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp &
      (cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(22, 4) = 0
H_A1(23, 4) = 0
H_A1(24, 4) = 0
H_A1(1, 5) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(2, 5) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(3, 5) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 5) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_A1(5, 5) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(6, 5) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(7, 5) = 0
H_A1(8, 5) = 0
H_A1(9, 5) = 0
H_A1(10, 5) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(11, 5) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(12, 5) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(13, 5) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(11))
H_A1(14, 5) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx &
      (0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(6)) + ( &
      1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(9)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1) &
      *exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1) &
      *exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1) &
      *exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8))
H_A1(15, 5) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(10))
H_A1(16, 5) = 0
H_A1(17, 5) = 0
H_A1(18, 5) = 0
H_A1(19, 5) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(20, 5) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(18)) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(21, 5) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(22, 5) = 0
H_A1(23, 5) = 0
H_A1(24, 5) = 0
H_A1(1, 6) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(2, 6) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(3, 6) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*k_y) - t_sigma*exp(-0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(4, 6) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(5, 6) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(6, 6) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(7, 6) = 0
H_A1(8, 6) = 0
H_A1(9, 6) = 0
H_A1(10, 6) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(11, 6) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(12, 6) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(13, 6) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      *conjg(xi(10)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg &
      (xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(6)) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(13)) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(14, 6) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(10))
H_A1(15, 6) = (1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1 &
      )*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0 &
      /2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0, &
      1)*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0 &
      /2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(6)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)* &
      k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(9))
H_A1(16, 6) = 0
H_A1(17, 6) = 0
H_A1(18, 6) = 0
H_A1(19, 6) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      cmplx(0,1)*k_y)*conjg(xi(1)) - exp(cmplx(0,1)*k_y)*conjg(xi(15)) &
      - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp &
      (cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(20, 6) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(21, 6) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(18)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(22, 6) = 0
H_A1(23, 6) = 0
H_A1(24, 6) = 0
H_A1(1, 7) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(2, 7) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(3, 7) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(4, 7) = 0
H_A1(5, 7) = 0
H_A1(6, 7) = 0
H_A1(7, 7) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(8, 7) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(9, 7) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(10, 7) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(11, 7) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(12, 7) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(13, 7) = 0
H_A1(14, 7) = 0
H_A1(15, 7) = 0
H_A1(16, 7) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*k_y)*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(18)) + exp(-cmplx(0,1)*k_y)*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(17, 7) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*k_y))*conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)* &
      k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4))
H_A1(18, 7) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*k_y)*conjg(xi(1)) + exp(-cmplx(0,1)*k_y)*conjg(xi(15 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(19, 7) = 0
H_A1(20, 7) = 0
H_A1(21, 7) = 0
H_A1(22, 7) = (1.0d0/2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*(cmplx(0,1)*( &
      -exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(9)) - cmplx(0, &
      1)*exp(-cmplx(0,1)*k_y)*conjg(xi(10)) - cmplx(0,1)*exp(-cmplx(0,1 &
      )*k_y)*conjg(xi(13)) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(7 &
      ))
H_A1(23, 7) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*((2.0d0 &
      /3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(7)) + (1.0d0/2.0d0)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(10))
H_A1(24, 7) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      *conjg(xi(11)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(10)) + (1.0d0/2.0d0)*(cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg &
      (xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(7)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(13)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(1, 8) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(2, 8) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(3, 8) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(4, 8) = 0
H_A1(5, 8) = 0
H_A1(6, 8) = 0
H_A1(7, 8) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(8, 8) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(9, 8) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_A1(10, 8) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(11, 8) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(12, 8) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(13, 8) = 0
H_A1(14, 8) = 0
H_A1(15, 8) = 0
H_A1(16, 8) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*k_y))*conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)* &
      k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4))
H_A1(17, 8) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(1)) + (1.0d0/2.0d0 &
      )*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp &
      (-cmplx(0,1)*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      *conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(4))
H_A1(18, 8) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + exp(-cmplx(0,1)*k_y))*conjg(xi(17)) + (1.0d0/2.0d0)*(exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      -cmplx(0,1)*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      *conjg(xi(2)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(19, 8) = 0
H_A1(20, 8) = 0
H_A1(21, 8) = 0
H_A1(22, 8) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*((2.0d0 &
      /3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(7)) + (1.0d0/2.0d0)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(11)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(10))
H_A1(23, 8) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(7)) + ( &
      1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (11)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (14)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi &
      (8)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg &
      (xi(9))
H_A1(24, 8) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg( &
      xi(7)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(14)) + (1.0d0/2.0d0 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(10)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(11))
H_A1(1, 9) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(2, 9) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(3, 9) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(4, 9) = 0
H_A1(5, 9) = 0
H_A1(6, 9) = 0
H_A1(7, 9) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_A1(8, 9) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(9, 9) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(10, 9) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(11, 9) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(12, 9) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*k_y) - t_sigma*exp(0.5d0 &
      *cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(13, 9) = 0
H_A1(14, 9) = 0
H_A1(15, 9) = 0
H_A1(16, 9) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*k_y)*conjg(xi(1)) + exp(-cmplx(0,1)*k_y)*conjg(xi(15 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(17, 9) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + exp(-cmplx(0,1)*k_y))*conjg(xi(17)) + (1.0d0/2.0d0)*(exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      -cmplx(0,1)*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      *conjg(xi(2)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))*conjg(xi(1)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(18, 9) = (1.0d0/2.0d0)*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))*conjg(xi(1)) + (1.0d0/2.0d0 &
      )*(exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      exp(-cmplx(0,1)*k_y))*conjg(xi(15)) + (1.0d0/2.0d0)*(exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)* &
      k_y))*conjg(xi(17)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))* &
      conjg(xi(5)) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))*conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(4))
H_A1(19, 9) = 0
H_A1(20, 9) = 0
H_A1(21, 9) = 0
H_A1(22, 9) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      *conjg(xi(11)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(10)) + (1.0d0/2.0d0)*(cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg &
      (xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(7)) - &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)*conjg(xi(13)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(23, 9) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg( &
      xi(7)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + (1.0d0/2.0d0)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(14)) + (1.0d0/2.0d0 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(10)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y))*conjg(xi(11))
H_A1(24, 9) = (1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(6)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))*conjg(xi(9)) + (1.0d0/2.0d0)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(11)) + (1.0d0/2.0d0 &
      )*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))* &
      conjg(xi(14)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(8))
H_A1(1, 10) = 0
H_A1(2, 10) = 0
H_A1(3, 10) = 0
H_A1(4, 10) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(5, 10) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(6, 10) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(7, 10) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 10) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(9, 10) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(10, 10) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(11, 10) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(12, 10) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(13, 10) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      cmplx(0,1)*k_y)*conjg(xi(16)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(18)) + exp(cmplx(0,1)*k_y)*conjg(xi(2)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(14, 10) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(15, 10) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      cmplx(0,1)*k_y)*conjg(xi(1)) - exp(cmplx(0,1)*k_y)*conjg(xi(15)) &
      - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(16, 10) = 0
H_A1(17, 10) = 0
H_A1(18, 10) = 0
H_A1(19, 10) = (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*(-cmplx(0,1)*( &
      -exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(9)) + cmplx(0, &
      1)*exp(cmplx(0,1)*k_y)*conjg(xi(10)) + cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(13)) + cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(7))
H_A1(20, 10) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(10))
H_A1(21, 10) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      *conjg(xi(11)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(10)) + (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg &
      (xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(13)) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(22, 10) = 0
H_A1(23, 10) = 0
H_A1(24, 10) = 0
H_A1(1, 11) = 0
H_A1(2, 11) = 0
H_A1(3, 11) = 0
H_A1(4, 11) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(5, 11) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(6, 11) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(7, 11) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(8, 11) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(9, 11) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(10, 11) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(11, 11) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(12, 11) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(13, 11) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(14, 11) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(18)) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(15, 11) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(16, 11) = 0
H_A1(17, 11) = 0
H_A1(18, 11) = 0
H_A1(19, 11) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(10))
H_A1(20, 11) = (1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1 &
      )*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8 &
      )) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(12)) + ( &
      1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(9))
H_A1(21, 11) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(11))
H_A1(22, 11) = 0
H_A1(23, 11) = 0
H_A1(24, 11) = 0
H_A1(1, 12) = 0
H_A1(2, 12) = 0
H_A1(3, 12) = 0
H_A1(4, 12) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(5, 12) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(6, 12) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(7, 12) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 12) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(9, 12) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*k_y) - t_sigma*exp(-0.5d0 &
      *cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(10, 12) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(11, 12) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(12, 12) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_A1(13, 12) = (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/2.0d0)*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      cmplx(0,1)*k_y)*conjg(xi(1)) - exp(cmplx(0,1)*k_y)*conjg(xi(15)) &
      - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(4))
H_A1(14, 12) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) + exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(15, 12) = (1.0d0/2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*(-exp(cmplx(0,1)*k_y) - exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*(exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(5)) - exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(18)) + exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4))
H_A1(16, 12) = 0
H_A1(17, 12) = 0
H_A1(18, 12) = 0
H_A1(19, 12) = (1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      *conjg(xi(11)) + (1.0d0/2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(10)) + (1.0d0/2.0d0)*(-cmplx(0,1)*(-exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*(-cmplx(0,1)*(-exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg &
      (xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg(xi(13)) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(8))
H_A1(20, 12) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(14)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) &
      + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(11))
H_A1(21, 12) = (1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(10)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(13)) + ( &
      1.0d0/2.0d0)*(-sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(7)) + ( &
      1.0d0/2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)))*conjg(xi(6)) + (1.0d0/2.0d0)*(-2.0d0/3.0d0* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(9)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)* &
      k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(11)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)* &
      k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(14)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)* &
      k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(8))
H_A1(22, 12) = 0
H_A1(23, 12) = 0
H_A1(24, 12) = 0
H_A1(1, 13) = 0
H_A1(2, 13) = 0
H_A1(3, 13) = 0
H_A1(4, 13) = -cmplx(0,1)*xi(10)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) - cmplx(0,1)*xi(7)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )))
H_A1(5, 13) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0 &
      )*xi(13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/ &
      2.0d0)*xi(7)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(8)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(6, 13) = (1.0d0/2.0d0)*xi(10)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + cmplx(0,1)*xi(8)*exp(-cmplx(0,1)*k_y)
H_A1(7, 13) = 0
H_A1(8, 13) = 0
H_A1(9, 13) = 0
H_A1(10, 13) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) - xi(16)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + xi(2)*exp(-cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(11, 13) = xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx &
      (0,1)*k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(12, 13) = xi(1)*exp(-cmplx(0,1)*k_y) - xi(15)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(4)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)* &
      (exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(13, 13) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(14, 13) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(15, 13) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx( &
      0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*B*cos(theta_B))
H_A1(16, 13) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(17, 13) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(18, 13) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(19, 13) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(20, 13) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(21, 13) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(22, 13) = 0
H_A1(23, 13) = 0
H_A1(24, 13) = 0
H_A1(1, 14) = 0
H_A1(2, 14) = 0
H_A1(3, 14) = 0
H_A1(4, 14) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0 &
      )*xi(13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/ &
      2.0d0)*xi(7)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(8)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(5, 14) = (1.0d0/2.0d0)*xi(10)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + ( &
      1.0d0/2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*(-sqrt(3.0d0) &
      *exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + ( &
      1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + ( &
      1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*(-sqrt(3.0d0)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0 &
      ,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + ( &
      1.0d0/2.0d0)*xi(8)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + ( &
      1.0d0/2.0d0)*xi(9)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))
H_A1(6, 14) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*( &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(8)*(-sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )))
H_A1(7, 14) = 0
H_A1(8, 14) = 0
H_A1(9, 14) = 0
H_A1(10, 14) = xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx &
      (0,1)*k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(11, 14) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(15)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - xi(16)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)* &
      (-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) + xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(12, 14) = xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(13, 14) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_A1(14, 14) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(15, 14) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(16, 14) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(17, 14) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(18, 14) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(19, 14) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(20, 14) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(21, 14) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(22, 14) = 0
H_A1(23, 14) = 0
H_A1(24, 14) = 0
H_A1(1, 15) = 0
H_A1(2, 15) = 0
H_A1(3, 15) = 0
H_A1(4, 15) = (1.0d0/2.0d0)*xi(10)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + cmplx(0,1)*xi(8)*exp(-cmplx(0,1)*k_y)
H_A1(5, 15) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*( &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(8)*(-sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )))
H_A1(6, 15) = (1.0d0/2.0d0)*xi(10)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*(sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*(sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(8)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(9)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))
H_A1(7, 15) = 0
H_A1(8, 15) = 0
H_A1(9, 15) = 0
H_A1(10, 15) = xi(1)*exp(-cmplx(0,1)*k_y) - xi(15)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(4)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)* &
      (exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(11, 15) = xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(12, 15) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(15)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - xi(16)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*( &
      -exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      -exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) + xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(13, 15) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(14, 15) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0, &
      1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_A1(15, 15) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(16, 15) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(17, 15) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(18, 15) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*k_y) + t_sigma*exp( &
      0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(19, 15) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(20, 15) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(21, 15) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(22, 15) = 0
H_A1(23, 15) = 0
H_A1(24, 15) = 0
H_A1(1, 16) = cmplx(0,1)*xi(10)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi( &
      11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + cmplx(0,1)*xi(7)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )))
H_A1(2, 16) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-2.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(3, 16) = (1.0d0/2.0d0)*xi(10)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) - cmplx(0,1)*xi(8)*exp(cmplx(0,1)*k_y)
H_A1(4, 16) = 0
H_A1(5, 16) = 0
H_A1(6, 16) = 0
H_A1(7, 16) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + xi(16)*exp &
      (cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + xi(2)*exp(cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(3)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 16) = xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(5)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(9, 16) = xi(1)*exp(cmplx(0,1)*k_y) + xi(15)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(10, 16) = 0
H_A1(11, 16) = 0
H_A1(12, 16) = 0
H_A1(13, 16) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(14, 16) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(15, 16) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(16, 16) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(17, 16) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(18, 16) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx( &
      0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*B*cos(theta_B))
H_A1(19, 16) = 0
H_A1(20, 16) = 0
H_A1(21, 16) = 0
H_A1(22, 16) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(23, 16) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(24, 16) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(1, 17) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-2.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(2, 17) = (1.0d0/2.0d0)*xi(10)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(12)* &
      ((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(13)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (14)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(7)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (8)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )))
H_A1(3, 17) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(7)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(8)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 17) = 0
H_A1(5, 17) = 0
H_A1(6, 17) = 0
H_A1(7, 17) = xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(5)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 17) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15)* &
      (exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0 &
      /2.0d0)*xi(18)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + xi(2)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3 &
      )*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(9, 17) = xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(10, 17) = 0
H_A1(11, 17) = 0
H_A1(12, 17) = 0
H_A1(13, 17) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(14, 17) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma* &
      (exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)))
H_A1(15, 17) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(16, 17) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_A1(17, 17) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(18, 17) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(19, 17) = 0
H_A1(20, 17) = 0
H_A1(21, 17) = 0
H_A1(22, 17) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(23, 17) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(24, 17) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(1, 18) = (1.0d0/2.0d0)*xi(10)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) - cmplx(0,1)*xi(8)*exp(cmplx(0,1)*k_y)
H_A1(2, 18) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(7)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(8)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(3, 18) = (1.0d0/2.0d0)*xi(10)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(12)* &
      (-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      14)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(7)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )))
H_A1(4, 18) = 0
H_A1(5, 18) = 0
H_A1(6, 18) = 0
H_A1(7, 18) = xi(1)*exp(cmplx(0,1)*k_y) + xi(15)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) &
      + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 18) = xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(9, 18) = (1.0d0/2.0d0)*xi(1)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15) &
      *(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + ( &
      1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + xi(2)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(3)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(4)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(10, 18) = 0
H_A1(11, 18) = 0
H_A1(12, 18) = 0
H_A1(13, 18) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(14, 18) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(15, 18) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*k_y) + t_sigma*exp( &
      -0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(16, 18) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(17, 18) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0, &
      1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_A1(18, 18) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(19, 18) = 0
H_A1(20, 18) = 0
H_A1(21, 18) = 0
H_A1(22, 18) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(23, 18) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(24, 18) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(1, 19) = 0
H_A1(2, 19) = 0
H_A1(3, 19) = 0
H_A1(4, 19) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) - xi(16)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - xi(2)*exp(-cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(5, 19) = -xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(-exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(6, 19) = -xi(1)*exp(-cmplx(0,1)*k_y) - xi(15)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) - 1.0d0/2.0d0*xi(4)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5 &
      )*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(7, 19) = 0
H_A1(8, 19) = 0
H_A1(9, 19) = 0
H_A1(10, 19) = -cmplx(0,1)*xi(10)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - cmplx(0,1)*xi(7)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*((1.0d0/3.0d0)*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(11, 19) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0 &
      /3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi &
      (13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi( &
      7)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + ( &
      1.0d0/2.0d0)*xi(8)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(12, 19) = (1.0d0/2.0d0)*xi(10)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + cmplx(0,1)*xi(8)* &
      exp(-cmplx(0,1)*k_y)
H_A1(13, 19) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(14, 19) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(15, 19) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(16, 19) = 0
H_A1(17, 19) = 0
H_A1(18, 19) = 0
H_A1(19, 19) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(20, 19) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_A1(21, 19) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(22, 19) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(23, 19) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(24, 19) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(1, 20) = 0
H_A1(2, 20) = 0
H_A1(3, 20) = 0
H_A1(4, 20) = -xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(-exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx &
      (0,1)*k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(5, 20) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      + (1.0d0/2.0d0)*xi(15)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - xi(16)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)* &
      (-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) &
      - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(6, 20) = -xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - &
      1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(7, 20) = 0
H_A1(8, 20) = 0
H_A1(9, 20) = 0
H_A1(10, 20) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0 &
      /3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi &
      (13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi( &
      7)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + ( &
      1.0d0/2.0d0)*xi(8)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(11, 20) = (1.0d0/2.0d0)*xi(10)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(12)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*(sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(14)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      )) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*( &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(8)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(9)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y))
H_A1(12, 20) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(8)*(sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )))
H_A1(13, 20) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(14, 20) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(15, 20) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(16, 20) = 0
H_A1(17, 20) = 0
H_A1(18, 20) = 0
H_A1(19, 20) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(20, 20) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(21, 20) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_A1(22, 20) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(23, 20) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(24, 20) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(1, 21) = 0
H_A1(2, 21) = 0
H_A1(3, 21) = 0
H_A1(4, 21) = -xi(1)*exp(-cmplx(0,1)*k_y) - xi(15)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/2.0d0*xi(18)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) - 1.0d0/2.0d0*xi(4)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5 &
      )*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(5, 21) = -xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(17)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(19)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) - exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - &
      1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(6, 21) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) + &
      (1.0d0/2.0d0)*xi(15)*(-exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - xi(16)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*( &
      -exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      -exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp( &
      -cmplx(0,1)*k_y)) - xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y)) - &
      1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(-cmplx(0,1)*k_y))
H_A1(7, 21) = 0
H_A1(8, 21) = 0
H_A1(9, 21) = 0
H_A1(10, 21) = (1.0d0/2.0d0)*xi(10)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) - cmplx(0,1)*xi(13)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + cmplx(0,1)*xi(8)* &
      exp(-cmplx(0,1)*k_y)
H_A1(11, 21) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) - cmplx(0,1)*exp(-cmplx(0,1)* &
      k_y)) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y)) + (1.0d0/2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0, &
      1)*k_y)) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (2.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*(( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(8)*(sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )))
H_A1(12, 21) = (1.0d0/2.0d0)*xi(10)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (11)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0 &
      /2.0d0)*xi(12)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(13)*(-sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(14)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0 &
      /2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(7)*(-sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(8)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0 &
      /2.0d0)*xi(9)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*k_y))
H_A1(13, 21) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(14, 21) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(15, 21) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(16, 21) = 0
H_A1(17, 21) = 0
H_A1(18, 21) = 0
H_A1(19, 21) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_A1(20, 21) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(21, 21) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(22, 21) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(23, 21) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_A1(24, 21) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*k_y) + t_sigma*exp( &
      0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(1, 22) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15)*( &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + xi(16)*exp &
      (cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - xi(2)*exp(cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(3)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(2, 22) = -xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx(0,1)*k_y) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(3, 22) = -xi(1)*exp(cmplx(0,1)*k_y) + xi(15)*exp(cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5) &
      *(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 22) = 0
H_A1(5, 22) = 0
H_A1(6, 22) = 0
H_A1(7, 22) = cmplx(0,1)*xi(10)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi( &
      11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + cmplx(0,1)*xi(7)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(8, 22) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((2.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp &
      (cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(9, 22) = (1.0d0/2.0d0)*xi(10)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) - cmplx(0,1)*xi(8)* &
      exp(cmplx(0,1)*k_y)
H_A1(10, 22) = 0
H_A1(11, 22) = 0
H_A1(12, 22) = 0
H_A1(13, 22) = 0
H_A1(14, 22) = 0
H_A1(15, 22) = 0
H_A1(16, 22) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(17, 22) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_A1(18, 22) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      (1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(19, 22) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(20, 22) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(21, 22) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(22, 22) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(23, 22) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_A1(24, 22) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(1, 23) = -xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx(0,1)*k_y) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(2, 23) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15) &
      *(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0 &
      /2.0d0)*xi(18)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) - xi(2)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3 &
      )*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*( &
      -exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)))
H_A1(3, 23) = -xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx(0,1)*k_y) - exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 23) = 0
H_A1(5, 23) = 0
H_A1(6, 23) = 0
H_A1(7, 23) = (1.0d0/2.0d0)*xi(10)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) - &
      2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((2.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)*exp &
      (cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(8, 23) = (1.0d0/2.0d0)*xi(10)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(12)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(13)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (14)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(7)*(-sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (8)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)*( &
      -2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      )))
H_A1(9, 23) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-2.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (7)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(10, 23) = 0
H_A1(11, 23) = 0
H_A1(12, 23) = 0
H_A1(13, 23) = 0
H_A1(14, 23) = 0
H_A1(15, 23) = 0
H_A1(16, 23) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_A1(17, 23) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(18, 23) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_A1(19, 23) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_A1(20, 23) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma* &
      (exp(cmplx(0,1)*k_y) + exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)))
H_A1(21, 23) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(22, 23) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(23, 23) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_A1(24, 23) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_A1(1, 24) = -xi(1)*exp(cmplx(0,1)*k_y) + xi(15)*exp(cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) + exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0 &
      *k_y))) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*( &
      -exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0 &
      *k_y)) - exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5) &
      *(-exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(2, 24) = -xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(16)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      17)*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*( &
      exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(2)*(-exp(cmplx(0,1)*k_y) - exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(3, 24) = (1.0d0/2.0d0)*xi(1)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(15 &
      )*(exp(cmplx(0,1)*k_y) + exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + ( &
      1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(19)*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) - xi(2)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(3)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/2.0d0*xi(4)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*(-exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_A1(4, 24) = 0
H_A1(5, 24) = 0
H_A1(6, 24) = 0
H_A1(7, 24) = (1.0d0/2.0d0)*xi(10)*(-2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (2.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 2.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(12)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + cmplx(0,1)*xi(13)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(14)*(-cmplx(0,1)*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx( &
      0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      )) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) - cmplx(0,1)*xi(8)* &
      exp(cmplx(0,1)*k_y)
H_A1(8, 24) = (1.0d0/2.0d0)*xi(10)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(11)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(12)*((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (13)*(sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(14)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y) + ( &
      2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)*(-2.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (7)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) - 2.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8)*( &
      -sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(9, 24) = (1.0d0/2.0d0)*xi(10)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      11)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(12)* &
      ((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(13)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      14)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(6)* &
      ((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))) + (1.0d0/2.0d0)*xi(7)*(sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(8 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp &
      (cmplx(0,1)*k_y) + (2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(9)* &
      ((2.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))
H_A1(10, 24) = 0
H_A1(11, 24) = 0
H_A1(12, 24) = 0
H_A1(13, 24) = 0
H_A1(14, 24) = 0
H_A1(15, 24) = 0
H_A1(16, 24) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(17, 24) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_A1(18, 24) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_A1(19, 24) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(20, 24) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_A1(21, 24) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*k_y) + t_sigma*exp( &
      -0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_A1(22, 24) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_A1(23, 24) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_A1(24, 24) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
