#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

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
    for(size_t i = 0; i < width * height; i++)
    {
        state.image_color[i] = make_pixel(0, 0, 0);
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
    data_vertex input{};
    auto output = new data_geometry[3];
    auto vertex_ptr = state.vertex_data;

    switch(type)
    {
        case render_type::triangle:
        {
            for(size_t i = 0, j = 0; i < state.num_vertices; i++, j++) 
            {
                output[j].data = vertex_ptr;
                input.data = vertex_ptr;
                state.vertex_shader(input, output[j], state.uniform_data);
                vertex_ptr += state.floats_per_vertex;
                if(j == 2)
                {
                    rasterize_triangle(state, (const data_geometry**) &output);
                    j = -1;
                }
            }

            break;
        }
        case render_type::indexed:
            break;
        case render_type::fan:
            break;
        case render_type::strip:
            break;
        default:
            break;
    }

    delete [] output;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    int x[3];
    int y[3];

    int x_min = 0;
    int x_max = 0;
    int y_min = 0;
    int y_max = 0;

    int half_width  = state.image_width  / 2.0f;
    int half_height = state.image_height / 2.0f;

    // Convert coordinates from in[] into NDC coordinates.
    for(int i = 0; i < 3; i++) 
    {
        int temp_x = static_cast<int>(half_width  * (*in)[i].gl_Position[0] + (half_width  - 0.5));
        int temp_y = static_cast<int>(half_height * (*in)[i].gl_Position[1] + (half_height - 0.5));

        x[i] = temp_x;
        y[i] = temp_y;
        
        // Draw vertices.
        state.image_color[temp_x + temp_y * state.image_width] = make_pixel(255, 255, 255);
    }

     // Fragment shading
    float* data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment fragment_data{data};
    data_output* frag_out = new data_output;

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
        x_max = state.image_width - 1;
    }
    if(y_min < 0)
    {
        y_min = 0;
    }
    if(y_max > state.image_height)
    {
        y_max = state.image_height - 1;
    }

    // If the barycentric weights of the pixel add to 1, draw it.
    for(int j = y_min; j < y_max + 1; j++)
    {
        for(int i = x_min; i < x_max + 1; i++) 
        {
            float alpha = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / area;
            float beta =  (0.5f * ((x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / area;
            float gamma = (0.5f * ((x[0]*y[1] - x[1]*y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / area;

            // If the pixel is inside the triangle...
            if (alpha >= 0 && beta >= 0 && gamma >= 0) 
            {

                for(int k = 0; k < state.floats_per_vertex; k++)
                {
                    switch(state.interp_rules[k])
                    {
                        case interp_type::flat:
                        {
                            fragment_data.data[k] = in[0]->data[k];
                            break;
                        }
                        case interp_type::smooth:
                        {
                
                            break;
                        }
                        case interp_type::noperspective:
                        {
                
                            break;
                        }
                        default:
                        break;      
                    }
                }
               
                state.fragment_shader(fragment_data, *frag_out, state.uniform_data);
                state.image_color[i + j * state.image_width] =
                                     make_pixel(static_cast<int>(frag_out->output_color[0] * 255),
                                                static_cast<int>(frag_out->output_color[1] * 255),
                                                static_cast<int>(frag_out->output_color[2] * 255));

            }
        }
    }

    delete [] data;
    delete frag_out;

}
