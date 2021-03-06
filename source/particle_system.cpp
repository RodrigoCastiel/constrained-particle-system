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

  // Initialize particle coordinates, velocities and mass matrix.
  // mLength = 1.5*M_PI*radius/N;
  mLength = totalLength/N;
  double arcTheta = totalLength/radius;
  double mass = 1.0/(N+1);

  mM = Eigen::MatrixXd::Zero(2*N+2, 2*N+2);
  for (int i = 0; i < 2*(N+1); i++)
  {
    mM(i, i) = mass;
  }

  // std::cout << "M = \n" << mM << std::endl;

  // Distribute particles along the input ring (circle).
  for (int i = 0; i < mNumParticles; i++)
  {
    double theta = arcTheta*static_cast<double>(i)/N;
    // mX(0, i) = radius * sin(theta) + Xc(0);
    // mX(1, i) = radius * cos(theta) + Xc(1);

    mX(0, i) = 0.0;
    mX(1, i) = -static_cast<double>(i)/N;

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
  static int frame = 0;

  if (!mPaused)
  {
    // std::cout << "FRAME ________________________ " << frame++ << std::endl;

    for (int i = 0; i < 20; i++)
    {
      ParticleSystem::ExplicitEulerIntegration();
    }

    dynamic_cast<ConnectorsMesh*>(mConnector->GetMesh())->Update(mX);
  }
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

void ParticleSystem::SelectParticle(const glm::vec3 & C, const glm::vec3 & ray, int state)
{
  if (state == GLUT_UP)
  {
    std::cout << "Releasing particle " << mSelectedParticle << std::endl;
    mSelectedParticle = -1;
    return;
  }

  float t = -C[2]/ray[2];
  glm::vec3 I = C + t*ray;

  for (int i = 0; i < mNumParticles; i++)
  {
    glm::vec3 P;
    P[0] = mX(0, i);
    P[1] = mX(1, i);
    P[2] = 0.0f;
    double dist = glm::length( I - P );

    if (dist <= mLength/2.0)
    {
      mSelectedParticle = i;
      const double k = 100;
      mMouseForce(0) = k * (I[0] - P[0]);
      mMouseForce(1) = k * (I[1] - P[1]);

      std::cout << "Selected Particle = " << mSelectedParticle << std::endl;
      std::cout << "Mouse force = \n" << mMouseForce << std::endl;

      return;
    }
  }

  mSelectedParticle = -1;
}

void ParticleSystem::DragParticle(const glm::vec3 & C, const glm::vec3 & ray)
{
  if (mSelectedParticle != -1)
  {
    float t = -C[2]/ray[2];
    glm::vec3 I = C + t*ray;

    glm::vec3 P;
    P[0] = mX(0, mSelectedParticle);
    P[1] = mX(1, mSelectedParticle);
    P[2] = 0.0f;

    const double k = 100;
    mMouseForce(0) = k * (I[0] - P[0]);
    mMouseForce(1) = k * (I[1] - P[1]);

    std::cout << "Dragging Particle = " << mSelectedParticle << std::endl;
    std::cout << "Mouse force = \n" << mMouseForce << std::endl;
  }
}

void ParticleSystem::ExplicitEulerIntegration()
{
  int N = mNumParticles - 1;
  Eigen::VectorXd C(N + 3);
  Eigen::MatrixXd Fext(2, N+1);
  Eigen::MatrixXd gradC(N+3, 2*N + 2);
  Eigen::MatrixXd grad_dC(N+3, 2*N + 2);

  // Compute C(x), dC/dx and d(C')/dx.
  ConstraintFunc(mX, C);
  ConstraintGrad(mX, gradC);
  DiffConstraintGrad(mX, mV, grad_dC);

  // Compute external forces.
  ExternalForces(mX, Fext);

  // Assemble and solve linear system to find x'' and lambda.
  // Aw = d.
  Eigen::MatrixXd A(3*N+5, 3*N+5);               //    A =  |      M      (dC/dx)^T  |
  A << mM, gradC.transpose(),                    //         |    dC/dx        0      |
       gradC, Eigen::MatrixXd::Zero(N+3, N+3);   //

  Eigen::VectorXd d(3*N + 5);
  Eigen::VectorXd v(2*N+2);    // v =  | vx |
  v << mV.row(0).transpose(),  //      | vy |
       mV.row(1).transpose();

  // Baumgarte Stabilization to avoid drift.
  double b = 5; // Damping parameter.  
  d << Fext.row(0).transpose(),                      // d =  |                    fx                     | 
       Fext.row(1).transpose(),                      //      |                    fy                     |
       -((grad_dC * v) + 2*b*(gradC * v) + b*b*C );  //      |  - ddC/dx * x' - 2b* dC/dx * x' - b^2 * C |

  Eigen::VectorXd w = A.colPivHouseholderQr().solve(d);
  Eigen::MatrixXd accel(2, N+1);
  Eigen::VectorXd lambda(N+3);

  accel << w.segment(0, N+1).transpose(),
           w.segment(N+1, N+1).transpose();

  lambda << w.segment(2*N+2, N+3);
  Eigen::VectorXd fc = (gradC.transpose() * lambda);  // Constraint force.
  Eigen::MatrixXd Fc(2, N+1);
  Fc << fc.segment(0, N+1).transpose(), 
        fc.segment(N+1, N+1).transpose();

  Eigen::VectorXd fext(2*N+2);

  fext << Fext.row(0).transpose(),
          Fext.row(1).transpose();

  // Integrate acceleration to find velocities and integrate velocities to update positions.
  double m = 1.0/(N+1);
  double dt = mTimeStep;
  mV += ( dt * accel );
  mX += ( dt * mV );

  // double maxFcX = Fc.row(0).maxCoeff();
  // double maxFcY = Fc.row(1).maxCoeff();
  // // std::cout << "Max(FcX) = " << maxFcX << std::endl;
  // // std::cout << "Max(FcY) = " << maxFcY << std::endl;  
  // // std::cout << "Max(C)       = " << C.maxCoeff() << std::endl;
  // // std::cout << "Max(gradC)   = " << gradC.maxCoeff() << std::endl;
  // // std::cout << "Max(grad_dC) = " << grad_dC.maxCoeff() << std::endl;

  // // if (std::max(maxFcX, maxFcY) > 1e5)
  // // {
  // //   std::cout << "WARNING Numerical instability.\n";
  // // }
}

void ParticleSystem::ExternalForces(const Eigen::MatrixXd& X, Eigen::MatrixXd& Fext)
{
  Fext.setZero();

  // Evaluate to every particle.
  for (int i = 0; i < mNumParticles; i++)
  {
    Fext.col(i) = Eigen::Vector2d(0.0, -1.0);
  }

  if (mSelectedParticle != -1)
  {
    Fext.col(mSelectedParticle) += mMouseForce;

  }

  //std::cout << "Fext = \n" << Fext << std::endl;
}

void ParticleSystem::ConstraintFunc(const Eigen::MatrixXd& X, Eigen::VectorXd& C)
{
  Eigen::VectorXd C_rigid, C_pin;
  double C_ring;

  ParticleSystem::ConstraintRigid(X, C_rigid);
  C_ring = ParticleSystem::ConstraintRing(X);
  ParticleSystem::ConstraintPin(X, C_pin);

  // C_rigid *= 0;
  // C_ring  *= 0;

  // Should always be 0.
  C << C_rigid,
       C_pin,
       C_ring;

  // std::cout << "| C(x) |^2 = \n" << C.squaredNorm() << std::endl;
}

void ParticleSystem::ConstraintGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC)
{
  int N = mNumParticles - 1;
  Eigen::VectorXd gradC_ring(2*N + 2);
  Eigen::MatrixXd gradC_pin(2, 2*N + 2);
  Eigen::MatrixXd gradC_rigid(N, 2*N+2);

  ParticleSystem::ConstraintRigidGrad(X, gradC_rigid);
  ParticleSystem::ConstraintRingGrad(X, gradC_ring);
  ParticleSystem::ConstraintPinGrad(X, gradC_pin);

  // gradC_rigid *= 0;
  // gradC_ring  *= 0;

  gradC << gradC_rigid,
           gradC_pin,
           gradC_ring.transpose();

  // std::cout << "gradC = \n" << gradC << std::endl;
}

void ParticleSystem::DiffConstraintGrad(const Eigen::MatrixXd& X, const Eigen::MatrixXd& V, 
                                        Eigen::MatrixXd& grad_dC)
{
  int N = mNumParticles - 1;
  Eigen::VectorXd grad_dC_ring(2*N + 2);
  Eigen::MatrixXd grad_dC_pin(2, 2*N + 2);
  Eigen::MatrixXd grad_dC_rigid(N, 2*N+2);

  ParticleSystem::DiffConstraintRigidGrad(X, V, grad_dC_rigid);
  ParticleSystem::DiffConstraintRingGrad( X, V, grad_dC_ring);
  ParticleSystem::DiffConstraintPinGrad(  X, V, grad_dC_pin);

  // grad_dC_rigid *= 0;
  // grad_dC_ring  *= 0;

  grad_dC << grad_dC_rigid,
             grad_dC_pin,
             grad_dC_ring.transpose();
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

void ParticleSystem::ConstraintRigidGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC_rigid)
{
  int N = mNumParticles - 1;
  gradC_rigid.setZero();

  for (int i = 0; i < N; i++)
  {
    // Main diagonal (with respect to x).     dC/dxi
    gradC_rigid(i, i)   = 2*( X(0, i) - X(0, i+1) );

    // Diagonal above main diagonal (with respect to x).    dC/dx(i+1)  
    gradC_rigid(i, i+1) = -gradC_rigid(i, i);  

    // Main diagonal (with respect to y).     dC/dyi
    gradC_rigid(i, i+N+1)   = 2*( X(1, i) - X(1, i+1) );

    // Diagonal above main diagonal (with respect to y).    dC/dy(i+1)  
    gradC_rigid(i, i+1+N+1) = -gradC_rigid(i, i+N+1);  
  }

  // std::cout << "gradC_rigid = \n" << gradC_rigid << std::endl;
}

void ParticleSystem::ConstraintPinGrad(const Eigen::MatrixXd& X, Eigen::MatrixXd& gradC_pin)
{
  int N = mNumParticles - 1;
  gradC_pin.setZero();
  gradC_pin(0, 0)   = 1;
  gradC_pin(1, N+1) = 1;
  // std::cout << "gradC_pin = \n" << gradC_pin << std::endl;
}

void ParticleSystem::ConstraintRingGrad(const Eigen::MatrixXd& X, Eigen::VectorXd& gradC_ring)
{
  int N = mNumParticles - 1;
  double dC_dxn = 2*(X(0, N) - mRingCenter(0));
  double dC_dyn = 2*(X(1, N) - mRingCenter(1));

  gradC_ring << Eigen::VectorXd::Zero(N), dC_dxn, Eigen::VectorXd::Zero(N), dC_dyn;
  // std::cout << "gradC_ring = \n" << gradC_ring.transpose() << std::endl;
}

void ParticleSystem::DiffConstraintRigidGrad(const Eigen::MatrixXd& X, const Eigen::MatrixXd& V, 
                                             Eigen::MatrixXd& grad_dC_rigid)
{
  int N = mNumParticles - 1;
  grad_dC_rigid.setZero();

  for (int i = 0; i < N; i++)
  {
    // Main diagonal (with respect to x).     dC/dxi
    grad_dC_rigid(i, i) = 2*( V(0, i) - V(0, i+1) );

    // Diagonal above main diagonal (with respect to x).    dC/dx(i+1)  
    grad_dC_rigid(i, i+1) = -2*( V(0, i) - V(0, i+1) );  

    // Main diagonal (with respect to y).     dC/dyi
    grad_dC_rigid(i, i+N+1) = 2*( V(1, i) - V(1, i+1) );

    // Diagonal above main diagonal (with respect to y).    dC/dy(i+1)  
    grad_dC_rigid(i, i+1+N+1) = -2*( V(1, i) - V(1, i+1) );  
  }

  // std::cout << "grad_dC_rigid = \n" << grad_dC_rigid << std::endl;
}

void ParticleSystem::DiffConstraintRingGrad( const Eigen::MatrixXd& X, const Eigen::MatrixXd& V, 
                                             Eigen::VectorXd& grad_dC_ring)
{
  int N = mNumParticles - 1;
  double ddC_dxn = 2*V(0, N);
  double ddC_dyn = 2*V(1, N);
  grad_dC_ring << Eigen::VectorXd::Zero(N), ddC_dxn, Eigen::VectorXd::Zero(N), ddC_dyn;
  // std::cout << "grad_dC_ring = \n" << grad_dC_ring.transpose() << std::endl;
}

void ParticleSystem::DiffConstraintPinGrad(  const Eigen::MatrixXd& X, const Eigen::MatrixXd& V, 
                                             Eigen::MatrixXd& grad_dC_pin)
{
  grad_dC_pin.setZero();
  // std::cout << "grad_dC_pin = \n" << grad_dC_pin << std::endl;
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