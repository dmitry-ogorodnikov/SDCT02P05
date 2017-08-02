#pragma once

namespace Constants {
  const size_t stateSize = 6;
  const size_t N = 10;
  const double dt = 0.1;
  const double Lf = 2.67;
  const double refV = 100 * 0.44704;
  const double latency = 0.1;

  //ref mphSI(100)
  const double wCte = 250;
  const double wEpsi = 250;
  const double wDCte = 1500;
  const double wDEpsi = 1500;
  const double wV = 0.1;
  const double wSteering = 1;
  const double wAcc = 1;
  const double wDsteering = 250;
  const double wDacc = 1;

  const size_t numVars = N * stateSize + (N - 1) * 2;
  const size_t numConstraints = N * stateSize;

  const size_t idxX = 0;
  const size_t idxY = idxX + N;
  const size_t idxPsi = idxY + N;
  const size_t idxV = idxPsi + N;
  const size_t idxCte = idxV + N;
  const size_t idxEpsi = idxCte + N;
  const size_t idxSteering = idxEpsi + N;
  const size_t idxAcc = idxSteering + N - 1;
};
