#version 150

struct Light
{
  vec3 position;  // Center coordinates.
  vec3 normal;    // Direction vector (currently not being used).

  float La;   // in [0, 1].
  float Ld;   // in [0, 1].
  float Ls;   // in [0, 1].
};

struct Material
{
  vec3 Ka;
  vec3 Kd;
  vec3 Ks;
};

in vec4 v_color;
in vec3 v_normal;
in vec2 v_tex_coord;
in vec4 f_pos;

out vec4 c;

// Uniforms.
uniform int tex_on;
uniform sampler2D tex;

uniform int light_on;
uniform Light light;
uniform Material material;

void main()
{
  // float d = min(length(f_pos)/100.0, 1);  // The depth is just the norm of f_pos (camera coordinates).
  // color = d*vec4(1.0, 1.0, 1.0, 1) + (1-d)*color;

  if (light_on == 1)
  {
    vec3 position = vec3(0, 3, 0);
    vec3 normal   = normalize(vec3(0, -1, 1));

    vec3 l = normalize(position - f_pos.xyz);  // Vector from frag to light.
    vec3 n = normalize(v_normal);                    //
    vec3 r = reflect(f_pos.xyz, n);                  // Compute reflection for camera ray.
    float alpha = 1.0f;
    // Note: since f_pos is in camera coordinates, the incident ray is simply f_pos.xyz.          

    float q = length(position - f_pos.xyz);
    float att = 1;///(1 + 1e-6*q*q);

    // Basic lighting - Phong model.
    vec3 Ia = material.Ka *(light.La);
    vec3 Id = material.Kd *(light.Ld * max(dot(l, n), 0.0))*att;
    vec3 Is = material.Ks *(light.Ls * pow(max(dot(l, r), 0.0), alpha))*att;

    //c = color * (Ia + Id + Is) * dot(light.normal, l);
    vec4 color = (1-tex_on)*vec4(1, 1, 1, 1) + (tex_on)*texture(tex, v_tex_coord);
    c = color * vec4((Ia + Id + Is), 1.0);
  }
  else
  {
    c = (1-tex_on)*v_color + (tex_on)*texture(tex, v_tex_coord);
  }
}