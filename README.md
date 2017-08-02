# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

##The Model

Used simulator data:
 * `px` - current vehicle x-position in global coordinate system;
 * `py` - current vehicle y-position in global coordinate system;
 * `psi` - current vehicle orientation in global coordinate system;
 * `v` - current vehicle velocity (I convert this value from mph to m/s);
 * `steering` - current vehicle steering angle;
 * `acc` - current vehicle throttle value;
 * `ptsx` - x-positions of waypoints in global coordinate system;
 * `ptsy` - y-positions of waypoints in global coordinate system;

State data:
 1. `px` - vehicle x-position in local (vehicle) coordinate system;
 2. `py` - vehicle y-position in local (vehicle) coordinate system;
 3. `psi` - vehicle orientation in local (vehicle) coordinate system;
 4. `v` - vehicle velocity;
 5. `cte` - cross track error (the error between the center of the road and the vehicle's position);
 6. `epsi` - orientation error (the desired orientation subtracted from the current orientation);

Actuators data:
 1. `steering` - vehicle steering angle, which is limited [-deg2rad(25), deg2rad(25)] ;
 2. `acc` - vehicle throttle value, which is limited [-1, 1];

Update equations:
 * x[t+1] = x[t] + v[t] * cos(psi[t]) * dt
 * y[t+1] = y[t] + v[t] * sin(psi[t]) * dt
 * psi[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
 * v[t+1] = v[t] + a[t] * dt
 * cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
 * epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

##The prediction horizon

I used the prediction horizon value is equal 1 sec. I believe that 1 second is enough for relatively high speeds. For example, at a vehicle speed of 50 mph, it passes about 22 meters in 1 second. I set `dt` (elapsed duration between timesteps) equal to the latency value, i.e. 100 milliseconds. When I set `dt` less than 0.1 sec, the control becomes unstable. The number of timesteps in the horizon `N` equal to 10.

##Polynomial Fitting

Before polynomial fitting I transform waypoints from the global coord system to the vehicle coord system.

`localWayptsx[i] = (globalWayptsx[i]-px)*cos(psi) + (globalWayptsy[i]-py)*sin(psi);`
`localWayptsy[i] = -(globalWayptsx[i]-px)*sin(psi) + (globalWayptsy[i]-py)*cos(psi);`

##Model Predictive Control

To account for the latency between actuations commands I predict the vehicle state, which will be at time of execution of commands (that is, I calculate the state, which will be in 0.1 second) by kinematic model. Then I send the predicted state to the MPC controller.

The trajectory cost is computed by using following parameters:
 * cte[i]^2 - to reduce the cross track error;
 * epsi[i]^2 - to reduce the orientation error;
 * (v[i] - refV)^2 - to limit maximum speed of vehicle;
 * (cte[i+1] - cte[i])^2 - to smooth the vehicle motion;
 * (epsi[i+1] - epsi[i])^2 - to smooth the vehicle motion;
 * (steering[i])^2 - to minimize the use of steering;
 * (acc[i])^2 - to minimize the use of throttle;
 * (steering[i+1] - steering[i])^2 - to minimize sharp turns;
 * (acc[i+1] - acc[i])^2 - to minimize sudden braking;

I also used weights for these parameters to regulate their influence on the cost of the trajectory.
