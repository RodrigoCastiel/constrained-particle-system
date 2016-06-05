/******************************************+
*                                          *
*  Author: Rodrigo Castiel                 *
*  Email: rcrs2@cin.ufpe.br                *
*         castielr@usc.edu                 *
*                                          *
+*******************************************/

#pragma once

#include <vector>
#include <functional>
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
  void Setup(int numParticles, const Eigen::Vector2d& Xc, const Eigen::Vector2d& Xpin,
             double radius = 0.5, double totalLength = 1.0);

  void Animate();
  void Render() const;

  // Physical differential equations + constraint computation.

  // Computes the summation of all forces exerted by force fields and mouse.
  void ExternalForces(const Eigen::MatrixXd& X, Eigen::MatrixXd& Fext);

  // Computes the constraints for the particle system. It has to be ZERO.
  // Input: x. Output: C(x) to C.
  void ConstraintFunc(const Eigen::MatrixXd& X, Eigen::VectorXd& C);

  // Computes the gradient of the constraint with respect to x.
  // Input: x. Output: dC/dx to gradC.
  void ConstraintGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC);


  void ConstraintRigid(const Eigen::MatrixXd& X, Eigen::VectorXd& C_rigid);
  void ConstraintPin(const Eigen::MatrixXd& X, Eigen::VectorXd& C_pin);
  double ConstraintRing(const Eigen::MatrixXd& X);

  void ConstraintRigidGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC_rigid);
  void ConstraintPinGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC_pin);
  void ConstraintRingGrad(const Eigen::MatrixXd& X, Eigen::VectorXd& gradC_ring);

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
  Eigen::Vector2d mPinPosition;
  double mRingRadius;
  double mLength;

};

class ConnectorsMesh : public Mesh
{
public:
  ConnectorsMesh() { }
  virtual ~ConnectorsMesh() { }

  void Load(const Eigen::Matrix2Xd& X, float r = 1, float g = 1, float b = 1);
  void Update(const Eigen::Matrix2Xd& X);
};
