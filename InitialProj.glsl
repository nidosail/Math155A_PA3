
// A vertex shader and a fragment shader

// ***************************
// Vertex Shader, named "vertexShader_P_MV"
//     Sets the position and color of a vertex.
//   The projection and modelview matrices are used to position the vertex.
//   It copies the color to "theColor" so that the fragment shader can access it.
// ***************************
#beginglsl vertexshader vertexShader_P_MV
#version 330 core
layout (location = 0) in vec3 vPos;	   // Position in attribute location 0
layout (location = 1) in vec3 vColor;  // Color in attribute location 1    
out vec3 theColor;                     // output a color to the fragment shader    
uniform mat4 projectionMatrix;         // The projection matrix    
uniform mat4 modelviewMatrix;          // The model-view matrix    
void main()    
{    
   gl_Position = projectionMatrix * modelviewMatrix * vec4(vPos.x, vPos.y, vPos.z, 1.0);    
   theColor = vColor;    
}
#endglsl 

// ***************************
// Fragment shader named "fragmentShader_simple"
//     All this does is pass through the color.
//     This is the simplest possible fragment shader.
// ***************************
#beginglsl fragmentshader fragmentShader_simple
// Set a general color using a fragment shader. (A "fragment" is a "pixel".)
//    The color value is passed in, obtained from the colors on the vertices.
//    Color values range from 0.0 to 1.0.
//    First three values are Red/Green/Blue (RGB).
//    Fourth color value (alpha) is 1.0, meaning there is no transparency.
#version 330 core    
in vec3 theColor;        // Color value came from the vertex shader (smoothed)     
out vec4 FragColor;      // Color that will be used for the fragment    
void main()    
{    
   FragColor = vec4(theColor, 1.0f);   // Add alpha value of 1.0.    
}
#endglsl

