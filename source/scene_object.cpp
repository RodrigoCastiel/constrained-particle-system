#include "scene_object.h"

void SceneObject::Render() const
{
  if (IsInitialized())
  {
    mPipelineProgram->SetModelMatrix(mModelMatrix);

    if (HasMaterial())
    {
      mMaterial->Bind(mProgramHandle);
    }

    GLuint texLoc = glGetUniformLocation(mProgramHandle, "tex_on");
    if (HasTexture())
    {
      mTexture->Bind(mProgramHandle);
      glUniform1i(texLoc, 1);
    }
    else
    {
      glUniform1i(texLoc, 0);
    }

    GLuint lightOnLoc = glGetUniformLocation(mProgramHandle, "light_on");
    glUniform1i(lightOnLoc, mUsingLighting);

    mMesh->Render();
  }
}

void SceneObject::Animate()
{
  // Update position and rotation from corresponding velocities.
  mPos += mLinVelocity;
  mRot += mRotVelocity;

  mModelMatrix.LoadIdentity();
  mModelMatrix.Translate(mPos[0], mPos[1], mPos[2]);
  mModelMatrix.Rotate(mRot[2], 0, 0, 1);
  mModelMatrix.Rotate(mRot[0], 1, 0, 0);
  mModelMatrix.Rotate(mRot[1], 0, 1, 0);
  mModelMatrix.Scale(mScale[0], mScale[1], mScale[2]);
}

SceneObject::~SceneObject()
{
  if (mIsMeshOwner)
  {
    delete mMesh;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void Material::Bind(GLuint programHandle)
{
  // Set Material.
  GLuint locKa = glGetUniformLocation(programHandle, "material.Ka");
  GLuint locKd = glGetUniformLocation(programHandle, "material.Kd");
  GLuint locKs = glGetUniformLocation(programHandle, "material.Ks");

  glUniform3f(locKa, mKa[0], mKa[1], mKa[2]);
  glUniform3f(locKd, mKd[0], mKd[1], mKd[2]);
  glUniform3f(locKs, mKs[0], mKs[1], mKs[2]);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void Texture::Load(int width, int height)
{
  if (mBuffer != 0)  // There was a texture previously loaded.
  {
    glDeleteTextures(1, &mBuffer);
  }

  glGenTextures(1, &mBuffer);
  glBindTexture(GL_TEXTURE_2D, mBuffer);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glTexImage2D( GL_TEXTURE_2D,  // Target.
                0,              // Detail level - original.
                GL_RGB,         // How the colors are stored.
                width,   // Width.
                height,  // Height.
                0,                    // Border must be 0. 
                GL_RGB,               // Format (?).  
                GL_UNSIGNED_BYTE,        // Data type.
                nullptr                  // Buffer address.
  );
}

void Texture::Load(ImageIO* source)
{
  if (mBuffer != 0)  // There was a texture previously loaded.
  {
    glDeleteTextures(1, &mBuffer);
  }

  glGenTextures(1, &mBuffer);
  glBindTexture(GL_TEXTURE_2D, mBuffer);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  int bytesPerPixel = source->getBytesPerPixel();

  glTexImage2D( GL_TEXTURE_2D,  // Target.
                0,              // Detail level - original.
                GL_RGB,  // How the colors are stored.
                source->getWidth(),   // Width.
                source->getHeight(),  // Height.
                0,                    // Border must be 0. 
                GL_RGB,  // Format (?).  
                GL_UNSIGNED_BYTE,               // Data type.
                source->getPixels()             // Buffer address.
  );
}

bool Texture::Load(const std::string& fileName)
{
  bool successful = false;
  ImageIO* source = new ImageIO();
  if (source->loadJPEG(fileName.c_str()) == ImageIO::OK)
  {
    Texture::Load(source);
    successful = true;
  }
  else
  {
    std::cerr << "WARNING Texture file in " << fileName << " could not be loaded.\n";
  }

  delete source;
  return successful;
}

