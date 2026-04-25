H_E(1, 1) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(2, 1) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(3, 1) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(4, 1) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(5, 1) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(6, 1) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(7, 1) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + &
      B*sin(theta_B)*cos(phi_B))
H_E(8, 1) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(9, 1) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(10, 1) = 0
H_E(11, 1) = 0
H_E(12, 1) = 0
H_E(13, 1) = 0
H_E(14, 1) = 0
H_E(15, 1) = 0
H_E(16, 1) = (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0, &
      1)*exp(-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx &
      (0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(31)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(34)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(32 &
      )) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(35)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(39)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(45)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(49 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(38)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(44)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(51 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx &
      (0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg &
      (xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      37)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(40 &
      )) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(43)) + &
      (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(50)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(53))
H_E(17, 1) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(40)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(31)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(41)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(32)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(45)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(44)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(43)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(47)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(52))
H_E(18, 1) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(28)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(29)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(35)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(40)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(53))
H_E(19, 1) = 0
H_E(20, 1) = 0
H_E(21, 1) = 0
H_E(22, 1) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(3)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(6)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(60)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(63)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(69)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(14)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      17)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(4)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(55)) - 1.0d0/2.0d0*exp(-cmplx &
      (0,1)*k_y)*conjg(xi(58)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg &
      (xi(61)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0 &
      /2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(71))
H_E(23, 1) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(66)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(11)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(60)) &
      - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(61)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1 &
      )*k_y)*conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg( &
      xi(7))
H_E(24, 1) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(6)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) + &
      (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(17)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(66)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(71))
H_E(1, 2) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_E(2, 2) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(3, 2) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(4, 2) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(5, 2) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(6, 2) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(7, 2) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(8, 2) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + &
      B*sin(theta_B)*cos(phi_B))
H_E(9, 2) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(10, 2) = 0
H_E(11, 2) = 0
H_E(12, 2) = 0
H_E(13, 2) = 0
H_E(14, 2) = 0
H_E(15, 2) = 0
H_E(16, 2) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(40)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(31)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(41)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(32)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(45)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(44)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(43)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(47)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(52))
H_E(17, 2) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp &
      (-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(37)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(40)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(43)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(38)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(41)) + ( &
      1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(44)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(31)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(33)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(50 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx &
      (0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg &
      (xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg( &
      xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi &
      (32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      35)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(46 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(49)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(52))
H_E(18, 2) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(43)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(34)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(44)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(35)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(31)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(39)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(33)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(38)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(37)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(54))
H_E(19, 2) = 0
H_E(20, 2) = 0
H_E(21, 2) = 0
H_E(22, 2) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(66)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(11)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(60)) &
      - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(61)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1 &
      )*k_y)*conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg( &
      xi(7))
H_E(23, 2) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(5)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(1)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(61)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(60)) &
      + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(63)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(9))
H_E(24, 2) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(13)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(8)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(69)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) + &
      (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(14)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(72)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(9))
H_E(1, 3) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(2, 3) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(3, 3) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(4, 3) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(5, 3) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(6, 3) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*k_y) - t_sigma*exp(0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(7, 3) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(8, 3) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(9, 3) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + &
      B*sin(theta_B)*cos(phi_B))
H_E(10, 3) = 0
H_E(11, 3) = 0
H_E(12, 3) = 0
H_E(13, 3) = 0
H_E(14, 3) = 0
H_E(15, 3) = 0
H_E(16, 3) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(28)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(29)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(35)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(40)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(53))
H_E(17, 3) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(43)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(34)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(44)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(35)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(31)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(39)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(33)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(38)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/2.0d0 &
      *exp(-cmplx(0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(37)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(54))
H_E(18, 3) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(27)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(48)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(54))
H_E(19, 3) = 0
H_E(20, 3) = 0
H_E(21, 3) = 0
H_E(22, 3) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(6)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) + &
      (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(17)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp(-cmplx(0,1) &
      *k_y)*conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(66)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(71))
H_E(23, 3) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(13)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(17)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(8)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(69)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) + &
      (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(14)) - 1.0d0/2.0d0* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(72)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(9))
H_E(24, 3) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(70)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx &
      (0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg &
      (xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(66)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(69)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(72))
H_E(1, 4) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(2, 4) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(3, 4) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(4, 4) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(5, 4) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(6, 4) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(7, 4) = 0
H_E(8, 4) = 0
H_E(9, 4) = 0
H_E(10, 4) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(11, 4) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(12, 4) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(13, 4) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0) &
      *(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y))* &
      conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(29)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(32)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(35)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(28)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(31)) + ( &
      1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(34)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(20)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      )*conjg(xi(37)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(40)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(53)) &
      + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(38)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(44)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(39)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(45)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(49 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52))
H_E(14, 4) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(32)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(41)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(31)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(40)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(20)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(47)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      )*conjg(xi(52)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(21)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(22)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(26)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(30)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(44)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(49 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(45)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54 &
      ))
H_E(15, 4) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(38)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(37)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(21)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(35)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(40)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx( &
      0,1)*k_y)*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1 &
      )*k_y)*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      )*conjg(xi(53)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(23)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(27)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(36)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(41)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(50 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52 &
      ))
H_E(16, 4) = 0
H_E(17, 4) = 0
H_E(18, 4) = 0
H_E(19, 4) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(14)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(17)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(4)) + &
      (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(55)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/2.0d0)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(61)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      65)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(7)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(71)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(62)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(66)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(72)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(13)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(3)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(57)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(60)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(63)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(70)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(9))
H_E(20, 4) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(11)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi( &
      60)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(61)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(69)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(7)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(9))
H_E(21, 4) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(17)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(4)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(58)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(66)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(71)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(1)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(72)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(60)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(70))
H_E(22, 4) = 0
H_E(23, 4) = 0
H_E(24, 4) = 0
H_E(1, 5) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(2, 5) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(3, 5) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(4, 5) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_E(5, 5) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(6, 5) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(7, 5) = 0
H_E(8, 5) = 0
H_E(9, 5) = 0
H_E(10, 5) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(11, 5) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(12, 5) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(13, 5) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(32)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(41)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(31)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(40)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(20)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(47)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      )*conjg(xi(52)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(21)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(22)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(26)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(30)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(44)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(49 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(45)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54 &
      ))
H_E(14, 5) = (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)* &
      sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y)) &
      *conjg(xi(45)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(41)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(44)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(40 &
      )) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(43)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(19)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      35)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(46 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(49)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(52)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(33)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(36)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(47)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(50 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(31)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54))
H_E(15, 5) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(35)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(44)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(34)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(43)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(19)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(37)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(50)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(54)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(20)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(24)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(25)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(33)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(38)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(31)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(39)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53 &
      ))
H_E(16, 5) = 0
H_E(17, 5) = 0
H_E(18, 5) = 0
H_E(19, 5) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(11)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi( &
      60)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(61)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(69)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(7)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(9))
H_E(20, 5) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(3)) - 1.0d0/2.0d0 &
      *exp(cmplx(0,1)*k_y)*conjg(xi(57)) - 1.0d0/2.0d0*exp(cmplx(0,1)* &
      k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(60 &
      )) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) - 1.0d0/2.0d0* &
      exp(cmplx(0,1)*k_y)*conjg(xi(9)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(58)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(61)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(7)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(56)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(8))
H_E(21, 5) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(14)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(55)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) + (1.0d0 &
      /2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0*exp(cmplx &
      (0,1)*k_y)*conjg(xi(72)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg( &
      xi(9)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(61)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(69)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(71)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8))
H_E(22, 5) = 0
H_E(23, 5) = 0
H_E(24, 5) = 0
H_E(1, 6) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)* &
      (-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(2, 6) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(3, 6) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*k_y) - t_sigma*exp(-0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(4, 6) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(5, 6) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(6, 6) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(7, 6) = 0
H_E(8, 6) = 0
H_E(9, 6) = 0
H_E(10, 6) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(11, 6) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(12, 6) = (1.0d0/4.0d0)*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(13, 6) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(38)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(37)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(21)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(35)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(40)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx( &
      0,1)*k_y)*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1 &
      )*k_y)*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      )*conjg(xi(53)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(23)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(27)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(36)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(41)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(50 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52 &
      ))
H_E(14, 6) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(35)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(44)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(34)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(43)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(19)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)* &
      exp(cmplx(0,1)*k_y)*conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(37)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(50)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(54)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(20)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(24)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(25)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(33)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(38)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(31)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(39)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53 &
      ))
H_E(15, 6) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(21)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(27)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx( &
      0,1)*k_y)*conjg(xi(48)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(51)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(54)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(22)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(25)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(46)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(52)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(20)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(23)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(26)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(47)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(50)) - 1.0d0 &
      /2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(53))
H_E(16, 6) = 0
H_E(17, 6) = 0
H_E(18, 6) = 0
H_E(19, 6) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(17)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(4)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(58)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(66)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(71)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(1)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(72)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(60)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(70))
H_E(20, 6) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(14)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(55)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) + (1.0d0 &
      /2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0*exp(cmplx &
      (0,1)*k_y)*conjg(xi(72)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg( &
      xi(9)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(61)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(69)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(71)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8))
H_E(21, 6) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(12)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(15)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi &
      (66)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(69)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(72)) - 1.0d0/2.0d0*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(13)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(11)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(17)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(71))
H_E(22, 6) = 0
H_E(23, 6) = 0
H_E(24, 6) = 0
H_E(1, 7) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(2, 7) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(3, 7) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(4, 7) = 0
H_E(5, 7) = 0
H_E(6, 7) = 0
H_E(7, 7) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(8, 7) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_E(9, 7) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(10, 7) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(11, 7) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(12, 7) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(13, 7) = 0
H_E(14, 7) = 0
H_E(15, 7) = 0
H_E(16, 7) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(3)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(6)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(60)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(63)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(69)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) - 1.0d0/2.0d0*exp(-cmplx( &
      0,1)*k_y)*conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg( &
      xi(14)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(17)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx( &
      0,1)*k_y)*conjg(xi(55)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg( &
      xi(58)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0/2.0d0*exp(-cmplx &
      (0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg &
      (xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(71))
H_E(17, 7) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(66)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      60)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(65)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(69)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(7))
H_E(18, 7) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(6)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) &
      - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(17)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(66)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(71))
H_E(19, 7) = 0
H_E(20, 7) = 0
H_E(21, 7) = 0
H_E(22, 7) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0 &
      )*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp &
      (-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0* &
      sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      ))*conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(28)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(31)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(34)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(29)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(32)) + ( &
      1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(35)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(39)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(45)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(49 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(38)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(44)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(37)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(40)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(53))
H_E(23, 7) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(31)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(40)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(32)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(41)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(45)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(44)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(47)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(52))
H_E(24, 7) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(37)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(38)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(35)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(40)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(53))
H_E(1, 8) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(2, 8) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(3, 8) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(4, 8) = 0
H_E(5, 8) = 0
H_E(6, 8) = 0
H_E(7, 8) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(8, 8) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(9, 8) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0, &
      1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_E(10, 8) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(11, 8) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(12, 8) = -1.0d0/2.0d0*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(13, 8) = 0
H_E(14, 8) = 0
H_E(15, 8) = 0
H_E(16, 8) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(9)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(12)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(66)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      60)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(65)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(69)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(7))
H_E(17, 8) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(5)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(59)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(62)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(8)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(1)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(4)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(61)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(3)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(9))
H_E(18, 8) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(13)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(17)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(8)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(14)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(72)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(9))
H_E(19, 8) = 0
H_E(20, 8) = 0
H_E(21, 8) = 0
H_E(22, 8) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(31)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(40)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(32)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(41)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(45)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(44)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(20)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(43)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(47)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(52))
H_E(23, 8) = (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)* &
      k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0, &
      1)*exp(-cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx &
      (0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(40)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(43)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(41 &
      )) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg(xi(44)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(31)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(51 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(33)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(47)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(50 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) &
      *conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(35)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(52))
H_E(24, 8) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(34)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(43)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(35)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(44)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(31)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(39)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(33)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(38)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(37)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(54))
H_E(1, 9) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(2, 9) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(3, 9) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(4, 9) = 0
H_E(5, 9) = 0
H_E(6, 9) = 0
H_E(7, 9) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx( &
      0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*B*cos(theta_B))
H_E(8, 9) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(9, 9) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(10, 9) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(11, 9) = (1.0d0/2.0d0)*t_rashba*(-exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(12, 9) = -1.0d0/2.0d0*t_pi*exp(-cmplx(0,1)*k_y) - t_sigma*exp(0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(13, 9) = 0
H_E(14, 9) = 0
H_E(15, 9) = 0
H_E(16, 9) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(2)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(6)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(1)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(72)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) &
      - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(17)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(4)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(66)) - 1.0d0/2.0d0*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(71))
H_E(17, 9) = -1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(13)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(17)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(3)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x - 0.5d0*k_y))*conjg(xi(57)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(62)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(8)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(2)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(56)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(61)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(7)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(1)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(14)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      63)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*k_y)*conjg(xi(68)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(72)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(9))
H_E(18, 9) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(11)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(14)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(17)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(65)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(68)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(71)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(10)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(16)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(64)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(67)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(70)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(12)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(18)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi( &
      66)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(69)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(72))
H_E(19, 9) = 0
H_E(20, 9) = 0
H_E(21, 9) = 0
H_E(22, 9) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(37)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(38)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(34)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(36)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(41)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(54 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(22)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(35)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(40)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(53))
H_E(23, 9) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y &
      ) + cmplx(0,1)*exp(-cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0 &
      )*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0 &
      *sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)))*conjg(xi(34)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)))*conjg(xi(43)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))* &
      conjg(xi(35)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))*conjg( &
      xi(44)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(21)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(31)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(39)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(24)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(33)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(38)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))*conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(19)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(23)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(37)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)*conjg(xi(54))
H_E(24, 9) = (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(20)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(23)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(47)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(50)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(53)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y))*conjg(xi(19)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(22)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))*conjg(xi(49)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))*conjg(xi(52)) + ( &
      1.0d0/2.0d0)*exp(-cmplx(0,1)*k_y)*conjg(xi(21)) + (1.0d0/2.0d0)* &
      exp(-cmplx(0,1)*k_y)*conjg(xi(24)) + (1.0d0/2.0d0)*exp(-cmplx(0,1 &
      )*k_y)*conjg(xi(27)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)* &
      k_y)*conjg(xi(48)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y &
      )*conjg(xi(51)) + (1.0d0/2.0d0)*cmplx(0,1)*exp(-cmplx(0,1)*k_y)* &
      conjg(xi(54))
H_E(1, 10) = 0
H_E(2, 10) = 0
H_E(3, 10) = 0
H_E(4, 10) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(5, 10) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(6, 10) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(7, 10) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(8, 10) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(9, 10) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(10, 10) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(11, 10) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(12, 10) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(13, 10) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(11)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(14)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi &
      (17)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(4)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(55)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(58)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(61)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(65)) &
      + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0* &
      exp(cmplx(0,1)*k_y)*conjg(xi(7)) + (1.0d0/2.0d0)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(71)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(12)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(15)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(18)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(5)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(66)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(69)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(72)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(8)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(10)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(6)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(60)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(9))
H_E(14, 10) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(60)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(61)) + &
      (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0/2.0d0*exp &
      (cmplx(0,1)*k_y)*conjg(xi(69)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(9))
H_E(15, 10) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(12)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(17)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(3)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg( &
      xi(4)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(58)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(66)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(71)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(1)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(72)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(60)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(70))
H_E(16, 10) = 0
H_E(17, 10) = 0
H_E(18, 10) = 0
H_E(19, 10) = (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)* &
      k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/ &
      2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1 &
      )*exp(cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1 &
      )*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx &
      (0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(29)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(32)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(35)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(28)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(31 &
      )) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(34)) - &
      1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/2.0d0*exp( &
      cmplx(0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(26)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(37)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      40)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(43)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(50)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(53)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(38)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(41)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(44)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(39)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(45)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(49 &
      )) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52))
H_E(20, 10) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(41)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(32)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(40)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(31)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(43)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(47)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(44)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(45)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54 &
      ))
H_E(21, 10) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(29)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(28)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(35)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(40)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(36)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(41)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52))
H_E(22, 10) = 0
H_E(23, 10) = 0
H_E(24, 10) = 0
H_E(1, 11) = 0
H_E(2, 11) = 0
H_E(3, 11) = 0
H_E(4, 11) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(5, 11) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(6, 11) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(7, 11) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(8, 11) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*t_sigma*(exp(cmplx(0,1)*k_y) + exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(9, 11) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(10, 11) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(11, 11) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(12, 11) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(13, 11) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(60)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(61)) + &
      (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(65)) - 1.0d0/2.0d0*exp &
      (cmplx(0,1)*k_y)*conjg(xi(69)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(7)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(12)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y))*conjg(xi(58)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)* &
      (0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(9))
H_E(14, 11) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(3)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(6)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg( &
      xi(60)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(9)) + (1.0d0/2.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(1 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(4)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(55)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(58)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(61)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(7)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(5)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(59)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8))
H_E(15, 11) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(14)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(55)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0*exp( &
      cmplx(0,1)*k_y)*conjg(xi(72)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(9)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(56)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(61)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(7)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(17)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(71)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8))
H_E(16, 11) = 0
H_E(17, 11) = 0
H_E(18, 11) = 0
H_E(19, 11) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(33)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(41)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(32)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(40)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(31)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(20)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(29)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(43)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(47)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(51)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(30)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(44)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(45)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - &
      1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(50)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54 &
      ))
H_E(20, 11) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0) &
      *(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(42)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y))* &
      conjg(xi(45)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(38)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(41)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1 &
      )*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(44)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(37)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1) &
      *(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(40)) + ( &
      1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1)*exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg(xi(43)) - &
      1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/2.0d0*exp( &
      cmplx(0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(29)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(32)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(35)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(46)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      49)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(52)) &
      - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(30)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(33)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(36)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(50 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(28)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(31)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(51 &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54))
H_E(21, 11) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(44)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(35)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(43)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(34)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(37)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(50)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(54)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(33)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(38)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(31)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(39)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53 &
      ))
H_E(22, 11) = 0
H_E(23, 11) = 0
H_E(24, 11) = 0
H_E(1, 12) = 0
H_E(2, 12) = 0
H_E(3, 12) = 0
H_E(4, 12) = (1.0d0/2.0d0)*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(5, 12) = (1.0d0/2.0d0)*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(6, 12) = (1.0d0/4.0d0)*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(7, 12) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(8, 12) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)))
H_E(9, 12) = -1.0d0/2.0d0*t_pi*exp(cmplx(0,1)*k_y) - t_sigma*exp(-0.5d0* &
      cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(10, 12) = (1.0d0/4.0d0)*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*lambda_soc + (1.0d0/2.0d0)*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(11, 12) = (1.0d0/4.0d0)*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*lambda_soc + (1.0d0/2.0d0)*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(12, 12) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) - 1.0d0/2.0d0*mu
H_E(13, 12) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(12)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(17)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(3)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg( &
      xi(4)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(57)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(58)) - 1.0d0/2.0d0*exp(cmplx( &
      0,1)*k_y)*conjg(xi(66)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(71)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(1)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(18)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(5)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(55)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(59)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(72)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(11)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(2)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(56)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(6)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(60)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(70))
H_E(14, 12) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(1)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(14)) + (1.0d0/2.0d0)*exp(cmplx &
      (0,1)*k_y)*conjg(xi(18)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(55)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(63)) + ( &
      1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(68)) - 1.0d0/2.0d0*exp( &
      cmplx(0,1)*k_y)*conjg(xi(72)) + (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(9)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(15)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(16)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(2)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(56)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(61)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(69)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(7)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(70)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(13)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(17)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(3)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(57)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(62)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(67)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(71)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(8))
H_E(15, 12) = (1.0d0/2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(12)) + (1.0d0/ &
      2.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(15)) + (1.0d0/2.0d0)*exp( &
      cmplx(0,1)*k_y)*conjg(xi(18)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)* &
      conjg(xi(66)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(69)) - &
      1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(72)) + (1.0d0/2.0d0)*exp &
      (-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(10 &
      )) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(13)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(16)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(64)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(67)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(70)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(11)) + (1.0d0/ &
      2.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(14)) + (1.0d0/2.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(17)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(65)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(68)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(71))
H_E(16, 12) = 0
H_E(17, 12) = 0
H_E(18, 12) = 0
H_E(19, 12) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(39)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(30)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(38)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(29)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(37)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(28)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(22)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(35)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(40)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg &
      (xi(48)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi( &
      49)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)*conjg(xi(53 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(36)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(41)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(50)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(54)) - &
      1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(34)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(42)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(51)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52))
H_E(20, 12) = (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y &
      ) - cmplx(0,1)*exp(cmplx(0,1)*k_y))*conjg(xi(45)) + (1.0d0/2.0d0) &
      *((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y))*conjg(xi(36)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt &
      (3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)))*conjg(xi(44)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)))*conjg(xi(35)) + (1.0d0/2.0d0)*(-1.0d0/3.0d0*sqrt(3.0d0)* &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx( &
      0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))* &
      conjg(xi(43)) + (1.0d0/2.0d0)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - cmplx(0,1) &
      *exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))*conjg( &
      xi(34)) - 1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(19)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(23)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(27)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1) &
      *k_y)*conjg(xi(32)) - 1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) &
      *conjg(xi(37)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(46)) + (1.0d0/6.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(50)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(54)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(24)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(33)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(38)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(51)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52 &
      )) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(21)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) + (1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(31)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(39)) - 1.0d0/6.0d0*sqrt(3.0d0)*exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(48)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y))*conjg(xi(49)) - 1.0d0/2.0d0*cmplx(0,1)*exp( &
      -cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53 &
      ))
H_E(21, 12) = -1.0d0/2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(21)) - 1.0d0/ &
      2.0d0*exp(cmplx(0,1)*k_y)*conjg(xi(24)) - 1.0d0/2.0d0*exp(cmplx(0 &
      ,1)*k_y)*conjg(xi(27)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)* &
      k_y)*conjg(xi(48)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)* &
      conjg(xi(51)) - 1.0d0/2.0d0*cmplx(0,1)*exp(cmplx(0,1)*k_y)*conjg( &
      xi(54)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y))*conjg(xi(19)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(22)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(25)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(46)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(49)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(52)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(20)) - 1.0d0/2.0d0*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(23)) - 1.0d0/ &
      2.0d0*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))* &
      conjg(xi(26)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(47)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))*conjg(xi(50)) - 1.0d0/2.0d0*cmplx(0,1)*exp(-cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y))*conjg(xi(53))
H_E(22, 12) = 0
H_E(23, 12) = 0
H_E(24, 12) = 0
H_E(1, 13) = 0
H_E(2, 13) = 0
H_E(3, 13) = 0
H_E(4, 13) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(26)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(31)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(32)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(33)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(34)*(-1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(36)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0 &
      ,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(37)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi &
      (38)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(39)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(40)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi &
      (41)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(42)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(43)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi &
      (44)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(46)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(50)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi( &
      51)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(53)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 13) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(24)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      43)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(44)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi( &
      48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      6.0d0)*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(6, 13) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(26)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt &
      (3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi( &
      36)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(40)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(42)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(47)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(48)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      6.0d0)*sqrt(3.0d0)*xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(53)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(7, 13) = 0
H_E(8, 13) = 0
H_E(9, 13) = 0
H_E(10, 13) = -1.0d0/2.0d0*xi(1)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      10)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(11)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(12)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(14)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(17)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      18)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(3)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(4)*exp( &
      -cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(58)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(6)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(61)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(62)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(63)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(64)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(65)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(66)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(68)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0) &
      *xi(69)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(7)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(71)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(72)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(8)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(9)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y))
H_E(11, 13) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(12)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(15)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)* &
      exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(61)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(62)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(63)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(66)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(67)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(68)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(69)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(12, 13) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(12)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(17)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(18)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(3)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp( &
      -cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(60)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(66)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(72)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(13, 13) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(14, 13) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(15, 13) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(16, 13) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(17, 13) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(18, 13) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(19, 13) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(20, 13) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(21, 13) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(22, 13) = 0
H_E(23, 13) = 0
H_E(24, 13) = 0
H_E(1, 14) = 0
H_E(2, 14) = 0
H_E(3, 14) = 0
H_E(4, 14) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(24)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      43)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(44)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi( &
      48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      6.0d0)*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 14) = (1.0d0/2.0d0)*xi(19)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi( &
      36)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt &
      (3.0d0)*xi(46)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)* &
      xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/6.0d0)*sqrt(3.0d0)*xi(48)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(6, 14) = (1.0d0/2.0d0)*xi(19)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y &
      ) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(37)*exp(-cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(39)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(50)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(53)*exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx &
      (0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(7, 14) = 0
H_E(8, 14) = 0
H_E(9, 14) = 0
H_E(10, 14) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(12)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(15)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(5)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)* &
      exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(61)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(62)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(63)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(66)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(67)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(68)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(69)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(11, 14) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(55)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)* &
      exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(8)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(9)*exp(-cmplx(0,1)*k_y)
H_E(12, 14) = -1.0d0/2.0d0*xi(1)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(14)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(17)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(18)* &
      exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(2)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(3)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(55)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(57)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(69)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(7)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(70)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(71)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(72)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)*exp(-cmplx(0,1)*k_y)
H_E(13, 14) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*B*cos(theta_B))
H_E(14, 14) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(15, 14) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(16, 14) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(17, 14) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(18, 14) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(19, 14) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(20, 14) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(21, 14) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(22, 14) = 0
H_E(23, 14) = 0
H_E(24, 14) = 0
H_E(1, 15) = 0
H_E(2, 15) = 0
H_E(3, 15) = 0
H_E(4, 15) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(-cmplx(0,1)*k_y) &
      + (1.0d0/2.0d0)*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(26)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt &
      (3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi( &
      36)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(40)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(42)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(47)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(48)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      6.0d0)*sqrt(3.0d0)*xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(53)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 15) = (1.0d0/2.0d0)*xi(19)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y &
      ) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(37)*exp(-cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(39)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(50)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(53)*exp(cmplx(0,1 &
      )*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx &
      (0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(6, 15) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(24)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0) &
      *cmplx(0,1)*xi(46)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0, &
      1)*xi(48)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*cmplx(0,1)*xi(49)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0, &
      1)*k_y) + (1.0d0/2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(7, 15) = 0
H_E(8, 15) = 0
H_E(9, 15) = 0
H_E(10, 15) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(12)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(17)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(18)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(3)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp( &
      -cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(60)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(66)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(72)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(11, 15) = -1.0d0/2.0d0*xi(1)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(14)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(17)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(18)* &
      exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(2)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(3)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(55)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(57)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(69)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(7)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(70)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(71)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(72)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)*exp(-cmplx(0,1)*k_y)
H_E(12, 15) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(12)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(15)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(17)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(18)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0* &
      xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(66)*exp(-cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(69)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(72)*exp &
      (-cmplx(0,1)*k_y)
H_E(13, 15) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(14, 15) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0,1 &
      )*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(15, 15) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(16, 15) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(17, 15) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(18, 15) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*k_y) + t_sigma*exp( &
      0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(19, 15) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(20, 15) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(21, 15) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(22, 15) = 0
H_E(23, 15) = 0
H_E(24, 15) = 0
H_E(1, 16) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      28)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (29)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      30)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(31)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(32)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(33)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(34)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(36)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(37)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(38 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*xi(39)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(40)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(42)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(cmplx(0,1) &
      *k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(45)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      51)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(53)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      54)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(2, 16) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(cmplx(0, &
      1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(45)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)* &
      xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(3, 16) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(26)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(28)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(29)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(30)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(34)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      36)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (38)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      39)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(40)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(42)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - 1.0d0/2.0d0*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(48)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)* &
      xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(4, 16) = 0
H_E(5, 16) = 0
H_E(6, 16) = 0
H_E(7, 16) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(10) &
      *exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(11)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(12)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(13 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /2.0d0*xi(14)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi( &
      16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(17)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(18)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(3)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*k_y) - 1.0d0 &
      /2.0d0*xi(5)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(55)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(57)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(58)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(59)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(6)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(60)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(61)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      66)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(69)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /2.0d0*xi(71)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(8 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(9)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y))
H_E(8, 16) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(12)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(15)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)* &
      xi(5)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(58)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)*exp &
      (cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(62)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(63)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(66)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(9, 16) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(12)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(17 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(2)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(3)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/2.0d0*xi(5)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(56)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(57)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi &
      (58)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(60)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(66)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(10, 16) = 0
H_E(11, 16) = 0
H_E(12, 16) = 0
H_E(13, 16) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(14, 16) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(15, 16) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(16, 16) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(17, 16) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(18, 16) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(19, 16) = 0
H_E(20, 16) = 0
H_E(21, 16) = 0
H_E(22, 16) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(23, 16) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(24, 16) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(1, 17) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(cmplx(0, &
      1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(45)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)* &
      xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(2, 17) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(20 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(34)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(36)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(37)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(38)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(39)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx &
      (0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(50)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(51)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - 1.0d0/6.0d0*sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(3, 17) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(20 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y &
      ) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(37)*exp(cmplx(0,1)* &
      k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(39)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      51)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(54)*exp(cmplx(0,1)*k_y)
H_E(4, 17) = 0
H_E(5, 17) = 0
H_E(6, 17) = 0
H_E(7, 17) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(12)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(15)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)* &
      xi(5)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(58)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)*exp &
      (cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(62)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(63)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(66)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(8, 17) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp &
      (cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(5)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(55)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(59)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(6)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/2.0d0)*xi(61)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(62)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(63)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)* &
      xi(8)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(9)*exp(cmplx(0,1)*k_y)
H_E(9, 17) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(13) &
      *exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(14)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(17)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(3)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(61)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(62)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(63)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi &
      (67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(69)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(7)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(72)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi &
      (9)*exp(cmplx(0,1)*k_y)
H_E(10, 17) = 0
H_E(11, 17) = 0
H_E(12, 17) = 0
H_E(13, 17) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(14, 17) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(15, 17) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(16, 17) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx( &
      0,1)*B*cos(theta_B))
H_E(17, 17) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(18, 17) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(19, 17) = 0
H_E(20, 17) = 0
H_E(21, 17) = 0
H_E(22, 17) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(23, 17) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(24, 17) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(1, 18) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(26)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(28)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(29)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(30)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(34)*exp(cmplx(0,1 &
      )*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      36)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (38)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      39)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(40)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(42)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - 1.0d0/2.0d0*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(48)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)* &
      xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(2, 18) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(20 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y &
      ) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(37)*exp(cmplx(0,1)* &
      k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(39)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0 &
      ,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0 &
      )*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      51)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(54)*exp(cmplx(0,1)*k_y)
H_E(3, 18) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(26 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /2.0d0*xi(27)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*cmplx(0,1)*xi(46) &
      *exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /2.0d0*cmplx(0,1)*xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(48)*exp(cmplx(0,1)* &
      k_y) - 1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0* &
      cmplx(0,1)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(53)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)* &
      xi(54)*exp(cmplx(0,1)*k_y)
H_E(4, 18) = 0
H_E(5, 18) = 0
H_E(6, 18) = 0
H_E(7, 18) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(12)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(17 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(2)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(3)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0, &
      1)*k_y) - 1.0d0/2.0d0*xi(5)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(56)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(57)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi &
      (58)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(60)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(66)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(8, 18) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(13) &
      *exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(14)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(15)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(17)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(3)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(61)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(62)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(63)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi &
      (67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(69)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(7)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(72)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi &
      (9)*exp(cmplx(0,1)*k_y)
H_E(9, 18) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(12)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(14)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(15)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(17)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(65)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(66)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(68)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(69)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(72)*exp(cmplx(0,1)*k_y)
H_E(10, 18) = 0
H_E(11, 18) = 0
H_E(12, 18) = 0
H_E(13, 18) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(14, 18) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(15, 18) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*k_y) + t_sigma*exp( &
      -0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(16, 18) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)* &
      cmplx(0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(17, 18) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0,1 &
      )*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(18, 18) = -1.0d0/4.0d0*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(19, 18) = 0
H_E(20, 18) = 0
H_E(21, 18) = 0
H_E(22, 18) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(23, 18) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(24, 18) = -1.0d0/4.0d0*g*mu_B*(-cmplx(0,1)*B*sin(phi_B)*sin(theta_B &
      ) + B*sin(theta_B)*cos(phi_B))
H_E(1, 19) = 0
H_E(2, 19) = 0
H_E(3, 19) = 0
H_E(4, 19) = (1.0d0/2.0d0)*xi(1)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (10)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(11)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(12)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(14)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0) &
      *xi(15)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(18)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(4)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(58)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(6)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(61)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(62)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(63)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(65)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0 &
      )*xi(66)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x &
      + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(69)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(7)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(70)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(71)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0) &
      *xi(72)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 19) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(12)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(61)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(62)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(63)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(65)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0 &
      )*xi(66)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(69)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(6, 19) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(12)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(17)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(18)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(3)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(60)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(66)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(72)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(7, 19) = 0
H_E(8, 19) = 0
H_E(9, 19) = 0
H_E(10, 19) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      28)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(31)*((1.0d0/3.0d0)*sqrt &
      (3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(32)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(33)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(34)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(36)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx( &
      0,1)*k_y) + cmplx(0,1)*exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(37)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi( &
      38)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/3.0d0*sqrt(3.0d0)*xi(39)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(40)*exp(-cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(42)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(45)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (51)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(53)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (54)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(11, 19) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(-cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx &
      (0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(45)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0,1)*k_y) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt &
      (3.0d0)*xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(12, 19) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      28)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(34)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*xi(36)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(37)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0 &
      /2.0d0)*xi(38)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (39)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(40)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(42)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(47)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(48)*exp(-cmplx(0,1)*k_y) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt &
      (3.0d0)*xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(53)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(13, 19) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(14, 19) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(15, 19) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(16, 19) = 0
H_E(17, 19) = 0
H_E(18, 19) = 0
H_E(19, 19) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(20, 19) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(21, 19) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_E(22, 19) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(23, 19) = -cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(24, 19) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(1, 20) = 0
H_E(2, 20) = 0
H_E(3, 20) = 0
H_E(4, 20) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(-cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(12)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(61)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(62)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(63)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(64)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0 &
      *k_y)) + (1.0d0/2.0d0)*xi(65)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0 &
      )*xi(66)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(67)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(68)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(69)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 20) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(2)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(3)*exp( &
      -cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(55)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(58)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      -cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(60)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(8)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(9)*exp(-cmplx(0,1)*k_y)
H_E(6, 20) = (1.0d0/2.0d0)*xi(1)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(14)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(15)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(17)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(18)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(72)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9) &
      *exp(-cmplx(0,1)*k_y)
H_E(7, 20) = 0
H_E(8, 20) = 0
H_E(9, 20) = 0
H_E(10, 20) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(-cmplx(0 &
      ,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx &
      (0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(43)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(44)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(45)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(47)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0,1)*k_y) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt &
      (3.0d0)*xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(11, 20) = -1.0d0/2.0d0*xi(19)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(-cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(27)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(28)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(29)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(34)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(36)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(37)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(38)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(39)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*(-1.0d0/3.0d0*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx &
      (0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + &
      (1.0d0/2.0d0)*xi(41)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (42)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + &
      cmplx(0,1)*exp(-cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46) &
      *exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(50)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(51)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/6.0d0*sqrt(3.0d0)*xi(52)*exp(-cmplx(0,1)*k_y) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1) &
      *(-0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(12, 20) = -1.0d0/2.0d0*xi(19)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(27)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(37)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(39)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (51)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(13, 20) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(14, 20) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(15, 20) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(16, 20) = 0
H_E(17, 20) = 0
H_E(18, 20) = 0
H_E(19, 20) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(20, 20) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(21, 20) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx( &
      0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_E(22, 20) = cmplx(0,1)*t_rashba*exp(0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(23, 20) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(24, 20) = (1.0d0/2.0d0)*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(1, 21) = 0
H_E(2, 21) = 0
H_E(3, 21) = 0
H_E(4, 21) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(12)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(17)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(18)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(3)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp( &
      -cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(56)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(57)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(60)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(66)*exp &
      (-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(70)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(71)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(72)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(5, 21) = (1.0d0/2.0d0)*xi(1)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (13)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(14)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(15)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(17)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(18)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)* &
      exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(61)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(62)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(63)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)*exp(-cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(72)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9) &
      *exp(-cmplx(0,1)*k_y)
H_E(6, 21) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(12)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(17 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(18)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(64)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(65)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(66)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(67)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(68)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(69)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      70)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(71)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(72)*exp(-cmplx(0,1)*k_y)
H_E(7, 21) = 0
H_E(8, 21) = 0
H_E(9, 21) = 0
H_E(10, 21) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(-cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(23)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp &
      (cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(26)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(27)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      28)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      29)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (30)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1) &
      *exp(-cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(34)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(35)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)* &
      sqrt(3.0d0)*xi(36)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(37)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0 &
      /2.0d0)*xi(38)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (39)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(40)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(42)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(47)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(48)*exp(-cmplx(0,1)*k_y) - 1.0d0 &
      /6.0d0*sqrt(3.0d0)*xi(49)*exp(-cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt &
      (3.0d0)*xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(53)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))
H_E(11, 21) = -1.0d0/2.0d0*xi(19)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi( &
      20)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(25)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(26)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      - 1.0d0/2.0d0*xi(27)*exp(-cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(31)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0* &
      k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(32)*exp(-cmplx(0,1)*k_y) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*(( &
      1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)* &
      ((1.0d0/3.0d0)*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)*exp( &
      -cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(37)*exp(-cmplx(0,1 &
      )*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0) &
      *xi(39)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*xi(43)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      44)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (45)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(-cmplx(0,1)*k_y) + cmplx(0,1)* &
      exp(-cmplx(0,1)*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(-cmplx &
      (0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(-cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi &
      (51)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(12, 21) = -1.0d0/2.0d0*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(20)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(21)*exp &
      (-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(23)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(24)*exp(-cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(25)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) - 1.0d0/2.0d0*xi(26 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(27)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(46)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) &
      + (1.0d0/2.0d0)*cmplx(0,1)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0, &
      1)*xi(48)*exp(-cmplx(0,1)*k_y) + (1.0d0/2.0d0)*cmplx(0,1)*xi(49)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*cmplx(0,1)*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1)*xi(51)*exp(-cmplx(0, &
      1)*k_y) + (1.0d0/2.0d0)*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*cmplx(0,1 &
      )*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y &
      )) + (1.0d0/2.0d0)*cmplx(0,1)*xi(54)*exp(-cmplx(0,1)*k_y)
H_E(13, 21) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(14, 21) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(15, 21) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(16, 21) = 0
H_E(17, 21) = 0
H_E(18, 21) = 0
H_E(19, 21) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0, &
      1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *B*cos(theta_B))
H_E(20, 21) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(21, 21) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(22, 21) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(23, 21) = -1.0d0/2.0d0*t_rashba*(-exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x + 0.5d0*k_y)) + exp(-cmplx(0,1)*k_y))
H_E(24, 21) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*k_y) + t_sigma*exp( &
      0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(1, 22) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi( &
      10)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(11)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(12)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(13)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(14)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(15)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(16)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(17)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(18)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(4)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx( &
      0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi( &
      55)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(58)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      6)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*xi(60)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      2.0d0*xi(62)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(63)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(66)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(67)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(68)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(69)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi( &
      7)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(8)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(9)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y))
H_E(2, 22) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(12)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(58)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(62)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(63)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(66)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(3, 22) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(12)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi &
      (17)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(2)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(3)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(56)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(57)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(58)*exp(cmplx( &
      0,1)*k_y) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(60)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(66)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(4, 22) = 0
H_E(5, 22) = 0
H_E(6, 22) = 0
H_E(7, 22) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(31)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(32)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(33)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(34)*(-1.0d0/ &
      3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)*(-1.0d0/3.0d0*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(36)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0, &
      1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(37)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      38)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(39)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(40)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      41)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(42)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(43)*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi( &
      44)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(46)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(48)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(50)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      6.0d0*sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1) &
      *k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(8, 22) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(24)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(cmplx(0,1) &
      *k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(43 &
      )*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(44)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(48 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0) &
      *sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(9, 22) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(36 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(40)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(42)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(48)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0) &
      *sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))
H_E(10, 22) = 0
H_E(11, 22) = 0
H_E(12, 22) = 0
H_E(13, 22) = 0
H_E(14, 22) = 0
H_E(15, 22) = 0
H_E(16, 22) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(17, 22) = -1.0d0/6.0d0*sqrt(6.0d0)*lambda_soc
H_E(18, 22) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(19, 22) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(0.86602540378443865d0* &
      k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(20, 22) = -cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(21, 22) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0, &
      1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(22, 22) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(23, 22) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/3.0d0*sqrt(6.0d0)*cmplx(0 &
      ,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*B*cos(theta_B))
H_E(24, 22) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx &
      (0,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx &
      (0,1)*B*cos(theta_B))
H_E(1, 23) = (1.0d0/2.0d0)*xi(10)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(11)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(12)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0, &
      1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(58)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(59)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(6)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(61)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(62)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(63)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(64)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(65)*exp(cmplx(0,1)*k_y) - &
      1.0d0/2.0d0*xi(66)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(67)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(68)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(69)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(7)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/2.0d0)*xi(8)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(9)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(2, 23) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(2)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(3)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(4)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(5)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(55)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(56)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(57)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(58)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(59)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(6)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(60)*exp(cmplx( &
      0,1)*k_y) + (1.0d0/2.0d0)*xi(61)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(62)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(63)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(7)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(8) &
      *exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(9)*exp(cmplx(0,1)*k_y)
H_E(3, 23) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi( &
      13)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(14)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(15)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(17)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(61)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(62)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(63)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(67)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(68)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(69)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(7)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(70)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(71)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(72)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*k_y)
H_E(4, 23) = 0
H_E(5, 23) = 0
H_E(6, 23) = 0
H_E(7, 23) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(24)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(cmplx(0,1) &
      *k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(31)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(32)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(33)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt(3.0d0 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(43 &
      )*exp(cmplx(0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(44)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      3.0d0)*sqrt(3.0d0)*xi(45)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(46)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt( &
      3.0d0)*xi(47)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(48 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0) &
      *sqrt(3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(8, 23) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (20)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(26)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(27)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(28)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(29)*exp(cmplx(0,1) &
      *k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(30)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0 &
      )*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) &
      - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y) - 1.0d0/ &
      3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx(0,1) &
      *(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(36 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(40)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(41)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(42)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - &
      cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/ &
      3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt( &
      3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      cmplx(0,1)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y &
      ))) + (1.0d0/2.0d0)*xi(45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0 &
      ,1)*k_y) - cmplx(0,1)*exp(cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(46)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      47)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(48)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(52)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi( &
      53)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/6.0d0)*sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(9, 23) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (20)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y) &
      - 1.0d0/3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(37)*exp(cmplx(0,1 &
      )*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(39)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(50)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(53)*exp(cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(54)*exp(cmplx(0,1)*k_y)
H_E(10, 23) = 0
H_E(11, 23) = 0
H_E(12, 23) = 0
H_E(13, 23) = 0
H_E(14, 23) = 0
H_E(15, 23) = 0
H_E(16, 23) = (1.0d0/6.0d0)*sqrt(6.0d0)*lambda_soc
H_E(17, 23) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(18, 23) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0),( &
      1.0d0/2.0d0)*sqrt(2.0d0)))
H_E(19, 23) = cmplx(0,1)*t_rashba*exp(-0.5d0*cmplx(0,1)*k_y)*sin( &
      0.86602540378443865d0*k_x)
H_E(20, 23) = (1.0d0/2.0d0)*t_pi*exp(-cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x + 0.5d0*k_y)) + (1.0d0/2.0d0)*t_sigma*(exp(cmplx(0,1)*k_y) + &
      exp(-cmplx(0,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(21, 23) = -1.0d0/2.0d0*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0, &
      1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(22, 23) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/3.0d0)*sqrt(6.0d0)* &
      cmplx(0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/3.0d0)*sqrt(3.0d0)* &
      cmplx(0,1)*B*cos(theta_B))
H_E(23, 23) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
H_E(24, 23) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) + (1.0d0/2.0d0)*sqrt(2.0d0)*cmplx( &
      0,1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0, &
      1)*B*cos(theta_B))
H_E(1, 24) = -1.0d0/2.0d0*xi(1)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(12)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(16)*exp(cmplx &
      (0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi &
      (17)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(18)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(2)*exp( &
      cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0 &
      *xi(3)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(4)*exp(cmplx(0,1)* &
      k_y) + (1.0d0/2.0d0)*xi(5)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(55)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(56)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(57)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(58)*exp(cmplx( &
      0,1)*k_y) - 1.0d0/2.0d0*xi(59)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(6)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(60)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(64)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(65)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(66)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(70)*exp(cmplx( &
      0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(71 &
      )*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(72)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))
H_E(2, 24) = (1.0d0/2.0d0)*xi(1)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi( &
      13)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(14)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(15)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(16)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*xi(17)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(18)*exp(cmplx(0,1)*k_y) + ( &
      1.0d0/2.0d0)*xi(2)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(3)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(55)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(56)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(57)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(61)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(62)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(63)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(67)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(68)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(69)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(7)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(70)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(71)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(72)* &
      exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(8)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(9)*exp( &
      cmplx(0,1)*k_y)
H_E(3, 24) = -1.0d0/2.0d0*xi(10)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(11)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(12)*exp( &
      cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(13)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(14)*exp &
      (cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*xi(15)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*xi(16)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*xi(17 &
      )*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0 &
      /2.0d0*xi(18)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(64)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(65)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(66)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)* &
      xi(67)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + &
      (1.0d0/2.0d0)*xi(68)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(69)*exp(cmplx(0,1)*k_y) + (1.0d0/ &
      2.0d0)*xi(70)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(71)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(72)* &
      exp(cmplx(0,1)*k_y)
H_E(4, 24) = 0
H_E(5, 24) = 0
H_E(6, 24) = 0
H_E(7, 24) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*k_y) + &
      (1.0d0/2.0d0)*xi(23)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(27)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(28)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (29)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      30)*(-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(34)*exp(cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/3.0d0*sqrt( &
      3.0d0)*xi(35)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(36 &
      )*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(37)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi &
      (38)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      39)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(40)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(41)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(42)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(46)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(47)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(48)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0) &
      *sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*k_y) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(50)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(51)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(52)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/6.0d0*sqrt(3.0d0)*xi(53)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0* &
      sqrt(3.0d0)*xi(54)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y))
H_E(8, 24) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi &
      (20)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + ( &
      1.0d0/2.0d0)*xi(21)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(24)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(25)* &
      exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0 &
      /2.0d0)*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(cmplx(0,1)*k_y) - 1.0d0/3.0d0* &
      sqrt(3.0d0)*xi(31)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/3.0d0*sqrt(3.0d0)*xi(32)*exp(cmplx(0,1)*k_y) &
      - 1.0d0/3.0d0*sqrt(3.0d0)*xi(33)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(34)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(35)* &
      (-1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*(0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi(36)*( &
      -1.0d0/3.0d0*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)*exp( &
      cmplx(0,1)*k_y)) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(37)*exp(cmplx(0,1 &
      )*k_y) + (1.0d0/3.0d0)*sqrt(3.0d0)*xi(38)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/3.0d0)*sqrt( &
      3.0d0)*xi(39)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(43)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx &
      (0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/ &
      2.0d0)*xi(44)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - cmplx(0,1)*exp(cmplx(0, &
      1)*(0.86602540378443865d0*k_x - 0.5d0*k_y))) + (1.0d0/2.0d0)*xi( &
      45)*((1.0d0/3.0d0)*sqrt(3.0d0)*exp(cmplx(0,1)*k_y) - cmplx(0,1)* &
      exp(cmplx(0,1)*k_y)) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(46)*exp(cmplx &
      (0,1)*k_y) + (1.0d0/6.0d0)*sqrt(3.0d0)*xi(47)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/6.0d0)*sqrt( &
      3.0d0)*xi(48)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(49)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/6.0d0*sqrt(3.0d0) &
      *xi(50)*exp(cmplx(0,1)*k_y) - 1.0d0/6.0d0*sqrt(3.0d0)*xi(51)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*(-0.86602540378443865d0* &
      k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(53)*exp(cmplx(0,1)* &
      (0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(54)*exp(cmplx(0,1)*k_y)
H_E(9, 24) = (1.0d0/2.0d0)*xi(19)*exp(cmplx(0,1)*(-0.86602540378443865d0 &
      *k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(20)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(21)* &
      exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(22)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/2.0d0)*xi(23)* &
      exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(24)*exp(cmplx(0,1)*k_y) + (1.0d0/2.0d0)*xi(25)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) + (1.0d0/ &
      2.0d0)*xi(26)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0* &
      k_y)) + (1.0d0/2.0d0)*xi(27)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0* &
      cmplx(0,1)*xi(46)*exp(cmplx(0,1)*(-0.86602540378443865d0*k_x - &
      0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(47)*exp(cmplx(0,1)*( &
      0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)* &
      xi(48)*exp(cmplx(0,1)*k_y) - 1.0d0/2.0d0*cmplx(0,1)*xi(49)*exp( &
      cmplx(0,1)*(-0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/ &
      2.0d0*cmplx(0,1)*xi(50)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x &
      - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1)*xi(51)*exp(cmplx(0,1)*k_y &
      ) - 1.0d0/2.0d0*cmplx(0,1)*xi(52)*exp(cmplx(0,1)*( &
      -0.86602540378443865d0*k_x - 0.5d0*k_y)) - 1.0d0/2.0d0*cmplx(0,1) &
      *xi(53)*exp(cmplx(0,1)*(0.86602540378443865d0*k_x - 0.5d0*k_y)) - &
      1.0d0/2.0d0*cmplx(0,1)*xi(54)*exp(cmplx(0,1)*k_y)
H_E(10, 24) = 0
H_E(11, 24) = 0
H_E(12, 24) = 0
H_E(13, 24) = 0
H_E(14, 24) = 0
H_E(15, 24) = 0
H_E(16, 24) = -1.0d0/2.0d0*lambda_soc*(cmplx(-1.0d0/6.0d0*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(17, 24) = -1.0d0/2.0d0*lambda_soc*(cmplx((1.0d0/6.0d0)*sqrt(6.0d0), &
      -1.0d0/2.0d0*sqrt(2.0d0)))
H_E(18, 24) = -1.0d0/4.0d0*g*mu_B*(cmplx(0,1)*B*sin(phi_B)*sin(theta_B) &
      + B*sin(theta_B)*cos(phi_B))
H_E(19, 24) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(20, 24) = (1.0d0/2.0d0)*t_rashba*(exp(cmplx(0,1)*k_y) - exp(-cmplx(0 &
      ,1)*(-0.86602540378443865d0*k_x + 0.5d0*k_y)))
H_E(21, 24) = (1.0d0/2.0d0)*t_pi*exp(cmplx(0,1)*k_y) + t_sigma*exp( &
      -0.5d0*cmplx(0,1)*k_y)*cos(0.86602540378443865d0*k_x)
H_E(22, 24) = -1.0d0/4.0d0*delta_tri - 1.0d0/6.0d0*sqrt(3.0d0)*cmplx(0,1 &
      )*lambda_soc - 1.0d0/2.0d0*mu_B*((1.0d0/6.0d0)*sqrt(6.0d0)*cmplx( &
      0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0, &
      1)*B*sin(theta_B)*cos(phi_B) - 1.0d0/3.0d0*sqrt(3.0d0)*cmplx(0,1) &
      *B*cos(theta_B))
H_E(23, 24) = -1.0d0/4.0d0*delta_tri + (1.0d0/6.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*lambda_soc - 1.0d0/2.0d0*mu_B*(-1.0d0/6.0d0*sqrt(6.0d0)*cmplx &
      (0,1)*B*sin(phi_B)*sin(theta_B) - 1.0d0/2.0d0*sqrt(2.0d0)*cmplx(0 &
      ,1)*B*sin(theta_B)*cos(phi_B) + (1.0d0/3.0d0)*sqrt(3.0d0)*cmplx(0 &
      ,1)*B*cos(theta_B))
H_E(24, 24) = (1.0d0/4.0d0)*B*g*mu_B*cos(theta_B) + (1.0d0/2.0d0)*mu
