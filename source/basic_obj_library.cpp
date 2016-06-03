#include "basic_obj_library.h"

///////////////////////////////////////////////////////////////////////////////////////////////////

void AxisMesh::Load()
{
  mNumVertices = 6;
  mNumIndices  = 6;
  mDrawMode = GL_LINES;

  GLfloat vertices[] = {0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f,
                        1.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f,
                        
                        0.0f, 0.0f, 0.0f,   0.0f, 1.0f, 0.0f,
                        0.0f, 1.0f, 0.0f,   0.0f, 1.0f, 0.0f,

                        0.0f, 0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
                        0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 1.0f
                       };
  GLuint indices[] = {0, 1, 2, 3, 4, 5};

  mVertexSize = 6;
  mHasColors  = true;
  mHasNormals = false;
  mHasTexCoord = false;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  memcpy(mVertices, vertices, sizeof(vertices));
  memcpy(mIndices,  indices,  sizeof(indices));

  mInitialized = true;
  Mesh::UploadGLBuffers();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void SquareMesh::Load()
{
  mNumVertices = 4;
  mNumIndices  = 4;
  mDrawMode = GL_TRIANGLE_STRIP;

  GLfloat vertices[] = {-1.0f, -1.0f, 0.0f,   0.0f, 0.0f,
                        -1.0f, +1.0f, 0.0f,   0.0f, 1.0f,
                        +1.0f, +1.0f, 0.0f,   1.0f, 1.0f,
                        +1.0f, -1.0f, 0.0f,   1.0f, 0.0f
                       };
  GLuint indices[] = {1, 0, 2, 3};

  mVertexSize = 5;
  mHasColors  = false;
  mHasNormals = false;
  mHasTexCoord = true;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  memcpy(mVertices, vertices, sizeof(vertices));
  memcpy(mIndices,  indices,  sizeof(indices));

  mInitialized = true;
  Mesh::UploadGLBuffers();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

#define INDEX(a, b) (((w) * (b)) + a)

void GridMesh::Load()
{
  int w = mWidth;
  int h = mHeight;

  mNumVertices = (w * h);
  mNumIndices = (2 * w * h);
  mDrawMode = GL_LINE_STRIP;

  mVertexSize = 6;
  mHasColors  = true;
  mHasNormals = false;
  mHasTexCoord = false;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  // Initialize vertices.
  for (int y = 0; y < h; y++)
  {
    for (int x = 0; x < w; x++)
    {
      GLfloat* position = PositionAt(INDEX(x, y));
      GLfloat* color = ColorAt(INDEX(x, y));
      position[0] = static_cast<float>(x - w/2)/w;
      position[1] = 0.0;
      position[2] = static_cast<float>(y - h/2)/h;
      color[0] = 0.45f;
      color[1] = 0.45f;
      color[2] = 0.45f;
    }
  }

  // Initialize indices.
  // Wire frame Element array - indices are written in a zig-zag pattern,
  // first horizontally and then vertically. It uses GL_LINE_STRIP. 
  int index = 0;
  int x = 0, y = 0;
  int dx = 1, dy = -1;

  for (y = 0; y < h; y++)  // Horizontally.
  {
    x = ((dx == 1) ? 0 : w-1);
    for (int k = 0; k < w; k++, x += dx)
      mIndices[index++] = INDEX(x, y);
    dx *= -1;
  }

  // Start from the last point to allow continuity in GL_LINE_STRIP.
  x = ((dx == 1) ? 0 : w-1);
  y = h-1;  
  dy = -1;
  for (int i = 0; i < w; i++, x += dx)  // Vertically.
  {
    y = ((dy == 1) ? 0 : h-1);
    for (int j = 0; j < h; j++, y += dy)
      mIndices[index++] = INDEX(x, y);
    dy *= -1;
  }

  mInitialized = true;
  Mesh::UploadGLBuffers();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void TexturedSphereMesh::Load(bool completeDome)
{
  int n = mDetailLevel + 1;
  int h = mDetailLevel + 1;

  mNumVertices = (n * h);
  mNumIndices  = 2*(h-2)*n + 2*n + 2*(h-2);
  mDrawMode = GL_TRIANGLE_STRIP;

  mVertexSize  = 8;
  mHasColors   = false;
  mHasNormals  = true;
  mHasTexCoord = true;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  // Initialize vertices.
  for (int v = 0; v < h; v++)
  {
    for (int u = 0; u < n; u++)
    {
      int index = (v*n + u);
      GLfloat theta_u = (2*M_PI * u) / (n-1);
      GLfloat theta_v = (M_PI * v) / (h-1);

      GLfloat* position = PositionAt(index);
      GLfloat* texCoord = TexCoordAt(index);
      GLfloat* normal = NormalAt(index);
      
      texCoord[0] = 1 - static_cast<float>(u) / (n-1);
      if (completeDome)
      {
        texCoord[1] = 1 - static_cast<float>(v) / (h-1);
      }
      else
      {
        texCoord[1] = 2 - 2*static_cast<float>(v) / (h-1);
      }
      
      position[0] = cos(theta_u) * sin(theta_v);
      position[2] = sin(theta_u) * sin(theta_v);
      position[1] = cos(theta_v);

      normal[0] = 2*position[0];
      normal[1] = 2*position[1];
      normal[2] = 2*position[2];
    }
  }

  // Initialize indices.
  int index = 0;
  for (int v = 0; v < h-1; v++)
  {
    // Zig-zag pattern: alternate between top and bottom.
    for (int u = 0; u < n; u++)
    {
      mIndices[index++] = ((v+0)*n + u);
      mIndices[index++] = ((v+1)*n + u);
    }

    // Triangle row transition: handle discontinuity.
    if (v < h-2)
    {
      // Repeat last vertex and the next row first vertex to generate 
      // two invalid triangles and get continuity in the mesh.
      mIndices[index++] = ((v+1)*n + (n-1)); //INDEX(this->width-1, y+1);
      mIndices[index++] = ((v+1)*n + 0);     //INDEX(0, y+1);
    }
  }

  mInitialized = true;
  Mesh::UploadGLBuffers();
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void TexturedTerrainMesh::Load(ImageIO* heightmap)
{
  int w = mWidth;
  int h = mHeight;

  mNumVertices = (w * h);
  mNumIndices = 2*(h-2)*w + 2*w + 2*(h-2);
  mDrawMode = GL_TRIANGLE_STRIP;

  mVertexSize = 8;
  mHasColors  = false;
  mHasNormals = true;
  mHasTexCoord = true;

  mVertices = new GLfloat [mVertexSize * mNumVertices];
  mIndices  = new GLuint  [mNumIndices];

  // Initialize vertices.
  for (int y = 0; y < h; y++)
  {
    for (int x = 0; x < w; x++)
    {
      GLfloat* position = PositionAt(INDEX(x, y));
      GLfloat* texCoord = TexCoordAt(INDEX(x, y));
      GLfloat* normal = NormalAt(INDEX(x, y));

      position[0] = static_cast<float>(x - w/2)/w;
      position[2] = static_cast<float>(y - h/2)/h;
      position[1] = (heightmap != nullptr) ? (heightmap->getPixel(x, y, 0)/4.0f - 255/4) : 0.0f;
  
      normal[0] = 0;
      normal[1] = 1;
      normal[2] = 0;

      texCoord[0] = static_cast<float>(x);
      texCoord[1] = static_cast<float>(y);
    }
  }

  // Initialize indices.
  int index = 0;
  for (int y = 0; y < h-1; y++)
  {
    // Zig-zag pattern: alternate between top and bottom.
    for (int x = 0; x < w; x++)
    {
      mIndices[index++] = INDEX(x, y);
      mIndices[index++] = INDEX(x, y+1);
    }

    // Triangle row transition: handle discontinuity.
    if (y < h-2)
    {
      // Repeat last vertex and the next row first vertex to generate 
      // two invalid triangles and get continuity in the mesh.
      mIndices[index++] = INDEX(w-1, y+1);
      mIndices[index++] = INDEX(0, y+1);
    }
  }

  mInitialized = true;
  Mesh::UploadGLBuffers();
}