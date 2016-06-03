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

class ConnectorsMesh;

class ParticleSystem
{
public:
  ParticleSystem(BasicPipelineProgram* pipelineProgram, GLuint programHandle) 
  : mPipelineProgram(pipelineProgram), mProgramHandle(programHandle)
  { }

  // Sets up the system according to the number of particles and circle position. 
  void Setup(int numParticles, const Eigen::Vector2d& Xc, double radius = 0.5);

  void Animate();
  void Render() const;

  ~ParticleSystem();

private:
  BasicPipelineProgram* mPipelineProgram;
  GLuint mProgramHandle;
  TexturedSphere* mParticleObj { nullptr };
  SceneObject* mConnector      { nullptr };
  SceneObject* mRing           { nullptr };

  // Physical Simulation data.
  int mNumParticles { 0 };
  Eigen::Matrix2Xd mX;  // Positions.
  Eigen::Matrix2Xd mV;  // Velocities.
  Eigen::Vector2d mRingCenter;
  double mRingRadius;
};

class ConnectorsMesh : public Mesh
{
public:
  ConnectorsMesh() { }
  virtual ~ConnectorsMesh() { }

  void Load(const Eigen::Matrix2Xd& X);
  void Update(const Eigen::Matrix2Xd& X);
};
