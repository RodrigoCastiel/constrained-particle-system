#include "particle_system.h"

using Eigen::Matrix2Xd;

void ParticleSystem::Setup(int numParticles, const Eigen::Vector2d& Xc, const Eigen::Vector2d& Xpin,
                           double radius, double totalLength)
{
  mNumParticles = ((numParticles > 1) ? numParticles : 11) + 1;
  int N = numParticles;

  mX = Matrix2Xd(2, mNumParticles);
  mV = Matrix2Xd(2, mNumParticles);
  double x_max = 0.5;

  // Initialize particle coordinates and velocities.
  // mLength = 1.5*M_PI*radius/N;
  mLength = totalLength/N;
  double arcTheta = totalLength/radius;

  // Distribute particles along the input ring (circle).
  for (int i = 0; i < mNumParticles; i++)
  {
    double theta = arcTheta*static_cast<double>(i)/N;
    mX(0, i) = radius * sin(theta) + Xc(0);
    mX(1, i) = radius * cos(theta) + Xc(1);

    std::cout << "x: " << mX(0, i) << ", y: " << mX(1, i) << std::endl;

    mV(0, i) = 0.0;
    mV(1, i) = 0.0;
  }

  // Load sphere to represent particle.
  mParticleObj = new TexturedSphere(mPipelineProgram, mProgramHandle);
  mParticleObj->SetScale(0.02f, 0.02f, 0.02f);
  mParticleObj->Load("");

  // Load connectors mesh.
  ConnectorsMesh* mesh = new ConnectorsMesh();
  mesh->SetProgramHandle(mProgramHandle);
  mesh->Load(mX, 0, 0, 1);

  mConnector = new SceneObject(mPipelineProgram, mProgramHandle);
  mConnector->SetMesh(mesh);
  mConnector->SetMeshOwner(true);

  // Initialize ring and load its mesh.
  int ns = 100;
  mRingCenter = Xc;
  mPinPosition = Xpin;
  mRingRadius = radius;
  Eigen::Matrix2Xd Xr(2, ns);
  for (int i = 0; i < ns; i++)
  {
    double theta = 2*M_PI*static_cast<double>(i)/(ns-1);
    Xr(0, i) = radius * cos(theta) + Xc(0);
    Xr(1, i) = radius * sin(theta) + Xc(1);
  }

  mesh = new ConnectorsMesh();
  mesh->SetProgramHandle(mProgramHandle);
  mesh->Load(Xr, 1, 1, 0);

  mRing = new SceneObject(mPipelineProgram, mProgramHandle);
  mRing->SetMesh(mesh);
  mRing->SetMeshOwner(true);
}

void ParticleSystem::Animate()
{
  int N = mNumParticles - 1;
  Eigen::VectorXd C(N + 3);
  Eigen::MatrixXd Fext(2, N+1);

  ConstraintFunc(mX, C);
  ExternalForces(mX, Fext);
}

void ParticleSystem::Render() const
{
  // Render one sphere for each particle.
  for (int i = 0; i < mNumParticles; i++)
  {
    mParticleObj->SetPosition(mX(0, i), mX(1, i), 0);
    mParticleObj->Animate();
    mParticleObj->Render();
    mConnector->Render();
    mRing->Render();
  }

}

void ParticleSystem::ExternalForces(const Eigen::MatrixXd& X, Eigen::MatrixXd& Fext)
{
  // Evaluate to every particle.
  for (int i = 0; i < mNumParticles; i++)
  {
    Fext.col(i) = Eigen::Vector2d(0.0, -1.0);
  }

  // TODO: add mouse force.

  std::cout << "Fext:\n" << Fext << std::endl;
}

void ParticleSystem::ConstraintFunc(const Eigen::MatrixXd& X, Eigen::VectorXd& C)
{
  Eigen::VectorXd C_rigid, C_pin;
  double C_ring;

  ParticleSystem::ConstraintRigid(mX, C_rigid);
  C_ring = ParticleSystem::ConstraintRing(mX);
  ParticleSystem::ConstraintPin(mX, C_pin);

  // std::cout << "C_rigid:\n" << C_rigid << std::endl;
  // std::cout << "C_ring:\n" << C_ring << std::endl;
  // std::cout << "C_pin:\n" << C_pin << std::endl;

  C << C_rigid,
       C_pin,
       C_ring;

  std::cout << "C(x):\n" << C << std::endl;
}

void ParticleSystem::ConstraintRigid(const Eigen::MatrixXd& X, Eigen::VectorXd& C_rigid)
{
  int N = mNumParticles-1;
  int H = X.rows();
  // C_rigid(x) = || X(i) - X(i+1) ||^2 - L^2;
  C_rigid = ((X.block(0, 0, H, N) - X.block(0, 1, H, N)).colwise().squaredNorm());
  C_rigid = C_rigid.array() - pow(mLength, 2.0);
}

double ParticleSystem::ConstraintRing (const Eigen::MatrixXd& X)
{
  int N = mNumParticles - 1;
  return  (X.col(N) - mRingCenter).squaredNorm() - pow(mRingRadius, 2.0);
}

void ParticleSystem::ConstraintPin(const Eigen::MatrixXd& X, Eigen::VectorXd& C_pin)
{
  C_pin = X.col(0) - mPinPosition;
}

ParticleSystem::~ParticleSystem()
{
  delete mRing;
  delete mConnector;
  delete mParticleObj;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void ConnectorsMesh::Load(const Eigen::Matrix2Xd& X, float r, float g, float b)
{
  int n = X.cols();
  mNumVertices = n;
  mNumIndices  = n;
  mDrawMode = GL_LINE_STRIP;

  mVertexSize = 6;
  mHasColors  = true;
  mHasNormals = false;
  mHasTexCoord = false;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  for (int i = 0; i < n; i++)
  {
    float* position = PositionAt(i);
    float* color    = ColorAt(i);
    position[0] = X(0, i);
    position[1] = X(1, i);
    position[2] = 0.0;

    color[0] = r;
    color[1] = g;
    color[2] = b;

    mIndices[i] = i;
  }

  mInitialized = true;
  Mesh::UploadGLBuffers();
}

void ConnectorsMesh::Update(const Eigen::Matrix2Xd& X)
{
  int n = mNumVertices;
  glm::vec3 x, dx;

  for (int i = 0; i < n; i++)
  {
    float* position = PositionAt(i);
    position[0] = X(0, i);
    position[1] = X(1, i);
    position[2] = 0.0;
  }

  // Update vertices to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, mVbo);
  glBufferData(GL_ARRAY_BUFFER, mVertexSize * mNumVertices * sizeof(GLfloat), mVertices, GL_STATIC_DRAW);
}