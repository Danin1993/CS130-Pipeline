#include <vector>
#include <limits>
#include "driver_state.h"


driver_state::driver_state() { }

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width  = width;
    state.image_height = height;
    state.image_color  = nullptr;
    state.image_depth  = nullptr;
    
    // Render the whole image black to start.
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];

    for(int i = 0; i < width * height; ++i)
    {
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = std::numeric_limits<float>::max();
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
//    const auto vertex = state.vertex_data;
    const data_geometry *output[3];
    data_geometry temp[3];
    data_vertex input[3];

    switch(type)
    {
        case render_type::triangle:
        {
            int triangles = state.num_vertices / 3;
            int k = 0;

            for(int i = 0; i < triangles; ++i) 
            {
                for(int j = 0; j < 3; ++j, k += state.floats_per_vertex)
                {
                    input[j].data = &state.vertex_data[k];
                    temp[j].data = input[j].data;
                    state.vertex_shader(input[j], temp[j], state.uniform_data);
                    output[j] = &temp[j];
                }
                
                clip_triangle(state, output, 0);
            }

            break;
        }
        case render_type::indexed:
        {
            int triangles = state.num_triangles * 3;
            for(int i = 0; i < triangles; i += 3)
            {
                for(int j = 0; j < 3; ++j)
                {
                    input[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
                    temp[j].data = input[j].data;
                    state.vertex_shader(input[j], temp[j], state.uniform_data);
                    output[j] = &temp[j];
                }

                    clip_triangle(state, output, 0);
            }

            break;
        }
        case render_type::fan:
        {
            for(int i = 0; i < state.num_vertices; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    if(j != 0)
                    {
                        input[j].data = state.vertex_data + ((state.floats_per_vertex) * (i + j));
                    }
                    else
                    {
                        input[j].data = state.vertex_data + (j * state.floats_per_vertex);
                    }

                    temp[j].data = input[j].data;
                    state.vertex_shader(input[j], temp[j], state.uniform_data);
                    output[j] = &temp[j];
                }

                clip_triangle(state, output, 0);
            }
            break;
        }
        case render_type::strip:
        {
            for(int i = 0; i < state.num_vertices - 2; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    input[j].data = &state.vertex_data[(i+j) * state.floats_per_vertex];
                    temp[j].data = input[j].data;
                    state.vertex_shader(input[j], temp[j], state.uniform_data);
                    output[j] = &temp[j];
                }

                clip_triangle(state, output, 0);
            }

            break;
        }
        default:
            break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
    if(face == 1)
    {
        rasterize_triangle(state, in);
        return;
    }
    else
    {
        const data_geometry *temp[3] = { in[0], in[1], in[2] };
        data_geometry first[3], second[3];
        float temp_A, temp_B1, temp_B2 = 0.0f;
        vec4 p1, p2;

        vec4 a = (*in)[0].gl_Position;
        vec4 b = (*in)[1].gl_Position;
        vec4 c = (*in)[2].gl_Position;

        // If all of the verts are on the screen, do nothing.
        if(a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3])
        {
            return;
        }
        else
        {
            if(a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3])
            {
                first[0].data = new float[state.floats_per_vertex];
                first[1] = *in[1];
                first[2] = *in[2];

                temp_B1 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
                temp_B2 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);

                p1 = temp_B1 * a + (1 - temp_B1) * b;
                p2 = temp_B2 * c + (1 - temp_B2) * a;

                for(int i = 0; i < state.floats_per_vertex; ++i)
                {
                    switch(state.interp_rules[i])
                    {
                        case interp_type::flat:
                        {
                            first[0].data[i] = (*in)[0].data[i];
                            break;
                        }
                        case interp_type::smooth:
                        {
                            first[0].data[i] = temp_B2 * (*in)[2].data[i] + (1 - temp_B2) * (*in)[0].data[i];
                            break;
                        }
                        case interp_type::noperspective:
                        {
                            temp_A = temp_B2 * (*in)[2].gl_Position[3] / (temp_B2 * (*in)[2].gl_Position[3] + 
                                               (1 - temp_B2) * (*in)[0].gl_Position[3]);
                            first[0].data[i] = temp_A * (*in)[2].data[i] + (1 - temp_A) * (*in)[0].data[i];
                            break;
                        }
                        default:
                            break;
                    }
                }
               
                first[0].gl_Position = p2;

                temp[0] = &first[0];
                temp[1] = &first[1];
                temp[2] = &first[2];

                clip_triangle(state, temp, face + 1);

                ////////////////////////////////////////////////////////////////////////////////////////////////


                second[0].data = new float[state.floats_per_vertex];
                second[1] = (*in)[1];
                second[2] = first[0];

                for(int i = 0; i < state.floats_per_vertex; ++i)
                {
                    switch(state.interp_rules[i])
                    {
                        case interp_type::flat:
                        {
                            second[0].data[i] = (*in)[0].data[i];
                            break;
                        }
                        case interp_type::smooth:
                        {
                            second[0].data[i] = temp_B1 * (*in)[0].data[i] + (1 - temp_B1) * (*in)[1].data[i];
                            break;
                        }
                        case interp_type::noperspective:
                        {
                            temp_A = temp_B1 * (*in)[0].gl_Position[3] / (temp_B1 * (*in)[0].gl_Position[3] + 
                                               (1 - temp_B1) * (*in)[1].gl_Position[3]);
                            second[0].data[i] = temp_A * (*in)[0].data[i] + (1 - temp_A) * (*in)[1].data[i];
                            break;
                        }
                        default:
                            break;
                    }
                }
               
                second[0].gl_Position = p2;

                temp[0] = &second[0];
                temp[1] = &second[1];
                temp[2] = &second[0];           
            }

            clip_triangle(state, temp, face + 1);
        }
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    float x[3];
    float y[3];
    float z[3];

    float x_min = 0;
    float x_max = 0;
    float y_min = 0;
    float y_max = 0;

    float half_width  = state.image_width  / 2.0f;
    float half_height = state.image_height / 2.0f;

    // Convert coordinates from in[] into NDC coordinates.
    for(int i = 0; i < 3; ++i) 
    {
        float temp_x = (half_width  * ((*in)[i].gl_Position[0]/(*in)[i].gl_Position[3]) + (half_width  - 0.5f));
        float temp_y = (half_height * ((*in)[i].gl_Position[1]/(*in)[i].gl_Position[3]) + (half_height - 0.5f));
        float temp_z = (half_width  * ((*in)[i].gl_Position[2]/(*in)[i].gl_Position[3]) + (half_width  - 0.5f));

        x[i] = temp_x;
        y[i] = temp_y;
        z[i] = temp_z;
    }

    // Fragment shading...
    float* data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment fragment_data{data};
    data_output frag_out;

    // Formula for the area of the triangle.
    float area = (0.5f * ((x[1]*y[2] - x[2]*y[1]) - (x[0]*y[2] - x[2]*y[0]) + (x[0]*y[1] - x[1]*y[0])));

    // Calculate the biggest x and y values in the triangle, and only visit those.  
    x_min = std::min(std::min(x[0], x[1]), x[2]);
    x_max = std::max(std::max(x[0], x[1]), x[2]);
    y_min = std::min(std::min(y[0], y[1]), y[2]);
    y_max = std::max(std::max(y[0], y[1]), y[2]);

    // Make sure the mins / maxs are still on the screen.
    if(x_min < 0)
    {
        x_min = 0;
    }
    if(x_max > state.image_width)
    {
        x_max = state.image_width;
    }
    if(y_min < 0)
    {
        y_min = 0;
    }
    if(y_max > state.image_height)
    {
        y_max = state.image_height;
    }

    // If the barycentric weights of the pixel add to 1, draw it.
    for(int j = y_min + 1; j < y_max; ++j)
    {
        for(int i = x_min + 1; i < x_max; ++i) 
        {
            float temp_alpha = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / area;
            float temp_beta =  (0.5f * ((x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / area;
            float temp_gamma = (0.5f * ((x[0]*y[1] - x[1]*y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / area;

            // If the pixel is inside the triangle...
            if (temp_alpha >= 0 && temp_beta >= 0 && temp_gamma >= 0) 
            { 
                float alpha = temp_alpha;
                float beta = temp_beta;
                float gamma = temp_gamma;

                float temp_z = alpha * z[0] + beta * z[1] + gamma * z[2];

                if(temp_z < state.image_depth[i + j * state.image_width])
                {
                    state.image_depth[i + j * state.image_width] = temp_z;

                    for(int k = 0; k < state.floats_per_vertex; k++)
                    {
                        float k_gour = 0;

                        switch(state.interp_rules[k])
                        {
                            case interp_type::flat:
                            {
                                fragment_data.data[k] = (*in)[0].data[k];
                                break;
                            }
                            case interp_type::smooth:
                            {
                                k_gour = (alpha / (*in)[0].gl_Position[3]) +
                                         (beta  / (*in)[1].gl_Position[3]) +
                                         (gamma / (*in)[2].gl_Position[3]);

                                temp_alpha = alpha / (k_gour * (*in)[0].gl_Position[3]);
                                temp_beta  = beta  / (k_gour * (*in)[1].gl_Position[3]);
                                temp_gamma = gamma / (k_gour * (*in)[2].gl_Position[3]);

                              //  break;
                            }
                            case interp_type::noperspective:
                            {
                                fragment_data.data[k] = temp_alpha * (*in)[0].data[k] + 
                                                        temp_beta  * (*in)[1].data[k] + 
                                                        temp_gamma * (*in)[2].data[k];
                                break;
                            }
                            default:
                                break;      
                        }
                    }
               
                    state.fragment_shader(fragment_data, frag_out, state.uniform_data);
                    state.image_color[i + j * state.image_width] =
                                     make_pixel(static_cast<int>(frag_out.output_color[0] * 255),
                                                static_cast<int>(frag_out.output_color[1] * 255),
                                                static_cast<int>(frag_out.output_color[2] * 255));

                }
            }
        }
    }

    delete [] data;
}
