/******************************************+
*                                          *
*  Author: Rodrigo Castiel                 *
*  Email: rcrs2@cin.ufpe.br                *
*         castielr@usc.edu                 *
*                                          *
+*******************************************/

#pragma once

#include <Eigen/Dense>
#include "basic_obj_library.h"

class ParticleSystem
{
public:
  ParticleSystem(BasicPipelineProgram* pipelineProgram, GLuint programHandle) 
  : mBasicPipelineProgram(pipelineProgram), mProgramHandle(programHandle)
  { }

  void Setup(int numParticles);

  void Animate();
  void Render() const;

  ~ParticleSystem();

private:
  BasicPipelineProgram* mBasicPipelineProgram;
  GLuint mProgramHandle;

  // Physical Simulation data.
  int mNumParticles;
  Eigen::Matrix2Xd mX;  // Positions.
  Eigen::Matrix2Xd mV;  // Velocities.
};