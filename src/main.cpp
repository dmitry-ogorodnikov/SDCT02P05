#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "Constants.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

namespace
{
  // For converting back and forth between radians and degrees.
  constexpr double pi() { return M_PI; }
  double deg2rad(const double x) { return x * pi() / 180; }

  // Checks if the SocketIO event has JSON data.
  // If there is data the JSON object in string format will be returned,
  // else the empty string "" will be returned.
  std::string hasData(const std::string& s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != std::string::npos) {
      return "";
    }
    else if (b1 != std::string::npos && b2 != std::string::npos) {
      return s.substr(b1, b2 - b1 + 2);
    }
    return "";
  }

  // Function to convert from mph to m/s
  double mph2SI(const double velocity) {
    return velocity * 0.44704;
  }

  // Evaluate a polynomial.
  double polyeval(const Eigen::VectorXd& coeffs, const double x) {
    double result = 0.0;
    const size_t amountCoeffs = coeffs.size();
    for (size_t i = 0; i < amountCoeffs; ++i) {
      result += coeffs[i] * pow(x, i);
    }
    return result;
  }

  Eigen::VectorXd polyfit(const Eigen::Ref<Eigen::VectorXd> xvals, const Eigen::Ref<Eigen::VectorXd> yvals,
    const size_t order) {
    const size_t amountPoints = xvals.size();

    assert(amountPoints == static_cast<size_t>(yvals.size()));
    assert(order >= 1 && order <= amountPoints - 1);

    Eigen::MatrixXd A(amountPoints, order + 1);

    for (size_t i = 0; i < amountPoints; ++i) {
      A(i, 0) = 1.0;
    }

    for (size_t j = 0; j < amountPoints; ++j) {
      for (size_t i = 0; i < order; ++i) {
        A(j, i + 1) = A(j, i) * xvals(j);
      }
    }

    const auto Q = A.householderQr();
    const auto result = Q.solve(yvals);
    return result;
  }

  void transformWaypoints(Eigen::Ref<Eigen::VectorXd> ptsx, Eigen::Ref<Eigen::VectorXd>& ptsy,
    const double px, const double py, const double psi) {

    assert(ptsx.size() == ptsy.size());

    const size_t amountPoints = ptsx.size();

    for (size_t i = 0; amountPoints > i; ++i) {
      const double deltaX = ptsx[i] - px;
      const double deltaY = ptsy[i] - py;

      ptsx[i] = deltaX*cos(psi) + deltaY*sin(psi);
      ptsy[i] = -deltaX*sin(psi) + deltaY*cos(psi);
    }
  }
}

int main() {
  uWS::Hub h;

  //Display the MPC predicted trajectory 
  std::vector<double> mpc_x_vals;
  std::vector<double> mpc_y_vals;

  h.onMessage([&](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    std::string sdata = std::string(data).substr(0, length);
    //std::cout << sdata << std::endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      std::string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          std::vector<double> ptsx = j[1]["ptsx"];
          std::vector<double> ptsy = j[1]["ptsy"];

          //A Eigen::vector expression mapping an existing array of data
          Eigen::Ref<Eigen::VectorXd> wayptsx = Eigen::Map<Eigen::VectorXd>(ptsx.data(), ptsx.size());
          Eigen::Ref<Eigen::VectorXd> wayptsy = Eigen::Map<Eigen::VectorXd>(ptsy.data(), ptsy.size());

          const double px = j[1]["x"];
          const double py = j[1]["y"];
          const double psi = j[1]["psi"];
          const double v = mph2SI(j[1]["speed"]);
          const double steering = j[1]["steering_angle"];
          const double acc = j[1]["throttle"];
          
          //Transformation of waypoints from the global coord system to the vehicle coord system.
          transformWaypoints(wayptsx, wayptsy, px, py, psi);
          const auto coeffs = polyfit(wayptsx, wayptsy, 3);

          //state - (x, y, psi, v, cte, epsi)
          Eigen::VectorXd state = Eigen::VectorXd::Zero(Constants::stateSize);

          /*
           * x = x0 + v0 * cos(psi) * dt
           * y = y0 + v0 * sin(psi) * dt
           * psi = psi0 + v0 / Lf * delta * dt
           * v = v0 + a * dt
           * cte = f(x0) - y0 + v0 * sin(epsi0) * dt
           * epsi = psi0 - psides0 + v0 * delta / Lf * dt
           */
          state(0) = v * Constants::latency;
          state(1) = 0.0;
          state(2) = -v * steering * Constants::latency / Constants::Lf;
          state(3) = v + acc * Constants::latency;
          state(4) = polyeval(coeffs, 0) + v * sin(-atan(coeffs[1])) * Constants::latency;
          state(5) = -atan(coeffs[1]) + state(2);

          const auto solution = MPC::solve(state, coeffs, mpc_x_vals, mpc_y_vals);

          json msgJson;
          
          msgJson["steering_angle"] = solution[0]/(deg2rad(25)*Constants::Lf);
          msgJson["throttle"] = solution[1];

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = ptsx;
          msgJson["next_y"] = ptsy;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      }
      else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
    size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    }
    else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
    char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  }
  else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
