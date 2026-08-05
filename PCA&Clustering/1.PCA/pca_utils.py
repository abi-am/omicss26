from __future__ import division
import numpy as np
from bokeh.core.properties import Instance, String
from sklearn.decomposition import PCA
from bokeh.io import output_notebook, push_notebook,show
from bokeh.layouts import row, column
from bokeh.models import Slider,Range1d,ColumnDataSource, LayoutDOM
from bokeh.plotting import figure, show
from bokeh.util.compiler import TypeScript
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.graph_objs as go
import plotly.offline as py
import numpy as np
from ipywidgets import interactive, HBox, VBox
import plotly.express as px


X = np.array([[-0.83934975, -0.21160323],
       [ 0.67508491,  0.25113527],
       [-0.05495253,  0.36339613],
       [-0.57524042,  0.24450324],
       [ 0.58468572,  0.95337657],
       [ 0.5663363 ,  0.07555096],
       [-0.50228538, -0.65749982],
       [-0.14075593,  0.02713815],
       [ 0.2587186 , -0.26890678],
       [ 0.02775847, -0.77709049]])


import numpy as np
def orthogonal_projection(p,n):
    """
    Given a normal vector to a plane, n and a point of space p, computes the orthogonal projection of p into the plane
    Input:
        p: numpy array
        n: numpy array
    Output:
        numpy array
    """
    n = np.array(n)
    p = np.array(p)
    lambda_val = np.dot(p,n)/np.dot(n,n)
    return p - lambda_val * n

def orthogonal_set_projection(P,n):
    """
    Given a normal vector to a plane, n and a set of points in space P, computes the orthogonal projection of each point
    Input:
        P: numpy array (or any iterable) of points in space
        n: numpy array
    Output:
        numpy array of arrays with the orthogonal projections
    """
    l = []
    for p in P:
        l.append(orthogonal_projection(p,n))
    return np.array(l)

def plot_line(plt,f, domain = [-1,1], **kwargs):
    f_x = [f(x) for x in domain]
    return plt.line(domain,f_x, **kwargs)



def get_plane_base(P,n):
    a,b = n
    n = n/np.linalg.norm(n)
    if b == 0:
        v = np.array([0,1])
    else:
        v = np.array([1,-a/b])
    v = v/np.linalg.norm(v)
    m = np.array([v,n])
    P_changed = []
    for p in P:
        P_changed.append(m@p)
    return np.array(P_changed)
        

def random_point_circle(center=(0,0),radius=1,n=1):
    r = radius * np.sqrt(np.random.rand(n))
    theta = np.random.rand(n) * 2 * np.pi
    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    return np.array([x,y]).T


def rotation_matrix(angle):
    radians = angle*np.pi/180
    return np.matrix([[np.cos(radians), -np.sin(radians)], [np.sin(radians), np.cos(radians)]]).round(4)

def line_between_points(ax,p1,p2,**kwargs):
    return ax.line(x=[p1[0],p2[0]],y=[p1[1],p2[1]],**kwargs)

## Fits the PCA with 2 dimensions
pca = PCA(n_components = 2).fit(X)
pca_line = np.array([k*pca.components_[0] for k in [-1.5,1.5]])
orthogonal_set = orthogonal_set_projection(X,(0,1))

def plot_widget_2(doc):
    main_plot = figure(width=500, height=500)
    main_plot.grid.visible = False
    main_plot.xaxis.visible = False
    main_plot.yaxis.visible = False
    main_plot.outline_line_color = None 
    main_plot.toolbar.logo = None
    main_plot.toolbar_location = None
    main_plot.x_range=Range1d(-1.5, 1.5)
    main_plot.y_range=Range1d(-1.5, 1.5)
    ## Main scatter plot (static)
    main_scatter = main_plot.scatter(X[:,0],X[:,1],marker = 'o', size = 4, color = "#C00000")

    ## PCA line (static)
    pca_line_plot = line_between_points(main_plot,pca_line[0],pca_line[1], color = '#333333', line_width = 1.5, legend_label = "PCA Line")

    ## Rotating line initial state:
    rotating_line_initial_state = [(-1.5,0), (1.5,0)]

    rotating_line = line_between_points(main_plot,
                                        rotating_line_initial_state[0], 
                                        rotating_line_initial_state[1], 
                                        color = '#0096FF',
                                        line_width = 1.5
                                        )

    projection_points = main_plot.scatter(orthogonal_set[:,0],
                                          orthogonal_set[:,1],
                                          marker = 'x', 
                                          size = 10, 
                                          color = "#C00000",
                                          legend_label = "Projected Points")

    lines = []*len(X)
    for p,o in zip(X,orthogonal_set):
        lines.append(line_between_points(main_plot,p,o,color='#FF9300', line_width = 1))

    rhs_plot = figure(width=500, height=500)
    rhs_plot.grid.visible = False
    rhs_plot_points = get_plane_base(orthogonal_set,(0,1)).round(2)
    line_between_points(rhs_plot,[-1,0],[1.15,0], color = "#0096FF",line_width = 1.5)
    rhs_scatter = rhs_plot.scatter(rhs_plot_points[:,0], rhs_plot_points[:,1], color = "#C00000", size = 5, marker = "x")
    rhs_plot.xaxis.visible = False
    rhs_plot.yaxis.visible = False
    rhs_plot.outline_line_color = None
    rhs_plot.toolbar.logo = None
    rhs_plot.toolbar_location = None
    rhs_plot.x_range=Range1d(-1.5, 1.5)

    slider = Slider(
        title="Adjust rotation angle",
        start=0,
        end=360,
        value= 0,
        step=1,

    )
    data = rotating_line.data_source.data
    p0 = np.array([data['x'][0],data['y'][0]])
    p1 = np.array([data['x'][1],data['y'][1]])
    def update(attr,old,new):
        ang = new
        # Rotate both points
        p0r = rotation_matrix(ang)@p0
        p1r = rotation_matrix(ang)@p1
        # This is the normal vector to the rotated line. 
        if abs(p0r[0] - p1r[0]) < 1e-10:
            n_line = (1,0)
        else:
            n_line = (-(p0r[1] - p1r[1])/(p0r[0] - p1r[0]), 1)
        # The code below just update every plot we just created above
        orthogonal_to_line = orthogonal_set_projection(X,n_line)
        rotating_line.data_source.data['x'] = [p0r[0],p1r[0]]
        rotating_line.data_source.data['y'] = [p0r[1],p1r[1]]
        projection_points.data_source.data['x'] = orthogonal_to_line[:,0]
        projection_points.data_source.data['y'] = orthogonal_to_line[:,1]
        
        projection_plot_1d = get_plane_base(orthogonal_to_line,n_line)
        rhs_scatter.data_source.data['x'] = projection_plot_1d[:,0]
        
        for line,p,o in zip(lines, X,orthogonal_to_line):
            line.data_source.data['x'] = [p[0],o[0]]
            line.data_source.data['y'] = [p[1],o[1]]
        

    slider.on_change('value',update)
    doc.add_root(row(main_plot,column(slider,rhs_plot)))
    
def plot_3d_2d_graphs(X):
    df = pd.DataFrame(X, columns = ['x1','x2'])
    df['x3'] = df['x1'] + df['x2']
    fig = make_subplots(rows=1, cols=2, specs = [[{"type":"scatter3d"}, {"type":"scatter"}]])


    fig.add_trace(go.Scatter3d(x= df['x1'], 
                      y = df['x2'], 
                      z = df['x3'], 
                      mode = 'markers'),
                      row = 1,
                      col = 1
                        ).update_traces(marker = dict(color = "#C00000", symbol = "x", size = 2),
                                        ).update_layout(scene = dict(xaxis = dict(range = [-1.5,1.5], showgrid=False),
                                                                     xaxis_title = 'x1',
                                                                     yaxis = dict(range = [-1.5,1.5], showgrid = False),
                                                                     yaxis_title = 'x2',
                                                                     zaxis = dict(range = [-1.5,1.5], showgrid = False),
                                                                     zaxis_title = 'x3'))

    fig.add_trace(go.Scatter(x = df.rename(columns = {'x1':'z1'})['z1'],
                             y = df.rename(columns = {'x2':'z2'})['z2'],
                             mode = "markers",
                             showlegend=False
                            ),
                    row = 1,
                    col = 2)
    return fig
    fig.show()
    
def plot_widget():
    main_plot = px.scatter(x = X[:,0], y = X[:,1])

    main_plot.data[0]['marker']['color'] = "#C00000"
    main_plot.data[0]['marker']['symbol'] = "x-thin-open"  
    main_plot.data[0]['marker']['size'] = 10

    
    pca_line_plot = line_between_points(px,pca_line[0],pca_line[1])

    pca_line_plot.data[0]['line']['color'] = "#333333"
    pca_line_plot.data[0]['showlegend'] = True
    pca_line_plot.data[0]['name'] = 'PCA Line'
    ## Rotating line initial state:
    rotating_line_initial_state = [(-2,0), (2,0)]



    rotating_line = line_between_points(px,
                                            rotating_line_initial_state[0], 
                                            rotating_line_initial_state[1]
                                            )



    rotating_line.data[0]['line']['color'] = "#0096FF"
    rotating_line.data[0]['showlegend'] = True
    rotating_line.data[0]['name'] = 'Rotating Line'

    projection_plot = px.scatter(x = orthogonal_set[:,0], y = orthogonal_set[:,1])
    projection_plot.data[0]['showlegend'] = True
    projection_plot.data[0]['name'] = 'Projection Points'

    projection_plot.data[0]['marker']['color'] = "#C00000"
    


    projection_plot_1d = get_plane_base(orthogonal_set,(0,1))
    rhs_scatter = px.scatter(x=projection_plot_1d[:,0],y=projection_plot_1d[:,1])

    rhs_scatter.data[0]['marker']['color'] = "#FF9300"
    final_data = main_plot.data + pca_line_plot.data + rotating_line.data + projection_plot.data
    rhs_initial_point = [-1.5,0]
    rhs_final_point = [1.5,0]
    rhs_line = line_between_points(px,rhs_initial_point,rhs_final_point)
    rhs_line.data[0]['line']['color'] = "#0096FF"
    
    lines = []*len(X)
    for p,o in zip(X,orthogonal_set):
        line = line_between_points(px,p,o)
        line.data[0]['line']['color'] = "#FF9300"
        lines.append(line)
        final_data = line.data + final_data


    p0 = np.array([rotating_line.data[0]['x'][0],rotating_line.data[0]['y'][0]])
    p1 = np.array([rotating_line.data[0]['x'][1],rotating_line.data[0]['y'][1]])
    def update_orthogonal_line(i,p,o):
        fig.data[i]['x'] = np.array([p[0],o[0]])
        fig.data[i]['y'] = np.array([p[1],o[1]])

    n_line = (1,0)

    def update(angle):
        ang = angle
        with fig.batch_update():
            p0r = rotation_matrix(ang)@p0
            p1r = rotation_matrix(ang)@p1
            # This is the normal vector to the rotated line. 
            if abs(p0r[0] - p1r[0]) < 1e-10:
                n_line = (1,0)
            else:
                n_line = (-(p0r[1] - p1r[1])/(p0r[0] - p1r[0]), 1)

        # The code below just update every plot we just created above
            orthogonal_to_line = orthogonal_set_projection(X,n_line)
            dispatches_args = enumerate(zip(X,orthogonal_to_line))
            for i,(o,p) in dispatches_args:
                fig.data[i]['x'] = np.array([p[0],o[0]])
                fig.data[i]['y'] = np.array([p[1],o[1]])
            fig.data[-2]['x'] = np.array([p0r[0],p1r[0]])
            fig.data[-2]['y'] = np.array([p0r[1],p1r[1]])
            fig.data[-1]['x'] = orthogonal_to_line[:,0]
            fig.data[-1]['y'] = orthogonal_to_line[:,1]
            projection_plot_1d = get_plane_base(orthogonal_to_line,n_line)
            rhs_fig.data[1]['x'] = np.array(projection_plot_1d[:,0])

    #     projection_plot_1d = get_plane_base(orthogonal_to_line,n_line)
    #     rhs_scatter.data_source.data['x'] = projection_plot_1d[:,0]






    ## The figures must exist before `interactive` is built: `interactive` calls `update`
    ## once at construction time, and that call reads `fig` / `rhs_fig` from this closure.
    fig = go.FigureWidget(data = final_data ).update_yaxes(scaleanchor = 'x', scaleratio= 1, range = [-1,1], visible=False).update_xaxes(range = [-1.5,1.5], visible=False)
    rhs_fig = go.FigureWidget(data = rhs_line.data + rhs_scatter.data).update_yaxes(scaleanchor = 'x', scaleratio= 1, range = [-1,1], showgrid=False, visible=False).update_xaxes(range = [-1.5,1.5], showgrid=False, visible=False)

    rhs_fig.update_layout(dict(width = 500, height = 400, plot_bgcolor = 'rgba(0,0,0,0)', title="PCA Projection"))

    fig.update_layout(dict(width = 500, height = 570, plot_bgcolor = 'rgba(0,0,0,0)'))

    freq_slider = interactive(update, angle=(0, 180, 1))

    vb = HBox((fig,(VBox(( freq_slider,rhs_fig)))))
    #vb.layout.align_items = 'center'
    return vb


def _rotated_line_and_projections(angle):
    """
    Given a rotation angle in degrees, returns the endpoints of the rotated line, the
    orthogonal projections of X onto it, and the 1-dimensional coordinates of those
    projections along the line.
    """
    p0 = np.array([-1.5, 0.0])
    p1 = np.array([ 1.5, 0.0])
    R = rotation_matrix(angle)
    p0r = np.asarray(R @ p0).ravel()
    p1r = np.asarray(R @ p1).ravel()

    ## Normal vector to the rotated line
    if abs(p0r[0] - p1r[0]) < 1e-10:
        n_line = (1, 0)
    else:
        n_line = (-(p0r[1] - p1r[1])/(p0r[0] - p1r[0]), 1)

    projections = orthogonal_set_projection(X, n_line)
    coords_1d = get_plane_base(projections, n_line)
    return p0r, p1r, projections, coords_1d


def plot_widget_mpl():
    """
    Matplotlib + ipywidgets version of `plot_widget`. Same picture, but drawn into an
    `Image` widget behind a plain ipywidgets slider, both of which render in VS Code --
    unlike plotly `FigureWidget` objects nested inside an `HBox`.

    The figure is built once and only the artists' data is updated on each slider move,
    and the resulting PNG is swapped into a single long-lived `Image` widget. Nothing in
    the output is ever torn down, so the picture updates in place without blanking.
    """
    from io import BytesIO
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from ipywidgets import IntSlider, Image, VBox, Layout

    dpi = 100
    width_in, height_in = 11, 5

    ## A bare Figure + Agg canvas, deliberately not a pyplot figure: this keeps the widget
    ## out of the inline backend's figure registry, so it is never drawn a second time as
    ## a normal cell output and never needs closing.
    fig = Figure(figsize = (width_in, height_in), dpi = dpi)
    FigureCanvasAgg(fig)
    ax_left, ax_right = fig.subplots(1, 2)

    p0r, p1r, projections, coords_1d = _rotated_line_and_projections(0)

    ## Left panel: the plane, the rotating line and the projections onto it.
    ## The artists created here are kept and mutated later instead of being recreated.
    connectors = [ax_left.plot([p[0],o[0]], [p[1],o[1]], color = "#FF9300",
                               linewidth = 1)[0]
                  for p,o in zip(X, projections)]
    ax_left.plot(pca_line[:,0], pca_line[:,1], color = "#333333",
                 linewidth = 1.5, label = "PCA line")
    rotating_line, = ax_left.plot([p0r[0],p1r[0]], [p0r[1],p1r[1]], color = "#0096FF",
                                  linewidth = 1.5, label = "Rotating line")
    ax_left.scatter(X[:,0], X[:,1], marker = 'o', color = "#C00000",
                    s = 35, label = "Original points")
    projected_points = ax_left.scatter(projections[:,0], projections[:,1], marker = 'x',
                                       color = "#C00000", s = 55,
                                       label = "Projected points")
    ax_left.set_xlim(-1.5, 1.5)
    ax_left.set_ylim(-1.5, 1.5)
    ax_left.set_aspect('equal')
    ax_left.axis('off')
    ax_left.legend(loc = 'upper left', fontsize = 8, frameon = False)
    left_title = ax_left.set_title("Projection onto a line rotated by 0°")

    ## Right panel: the same points, now living on a single dimension
    ax_right.plot([-1.5,1.5], [0,0], color = "#0096FF", linewidth = 1.5)
    projected_1d = ax_right.scatter(coords_1d[:,0], np.zeros(len(coords_1d)),
                                    marker = 'x', color = "#C00000", s = 55)
    ax_right.set_xlim(-1.5, 1.5)
    ax_right.set_ylim(-1, 1)
    ax_right.set_aspect('equal')
    ax_right.axis('off')
    ax_right.set_title("PCA projection (1 dimension)")

    def png_bytes(angle):
        p0r, p1r, projections, coords_1d = _rotated_line_and_projections(angle)

        for line,p,o in zip(connectors, X, projections):
            line.set_data([p[0],o[0]], [p[1],o[1]])
        rotating_line.set_data([p0r[0],p1r[0]], [p0r[1],p1r[1]])
        projected_points.set_offsets(projections)
        projected_1d.set_offsets(np.column_stack([coords_1d[:,0],
                                                  np.zeros(len(coords_1d))]))
        left_title.set_text(f"Projection onto a line rotated by {angle}°")

        ## Light PNG compression: ~10 ms/frame cheaper to encode for ~9 KB more per
        ## frame, which is the right trade when the kernel is local.
        buffer = BytesIO()
        fig.savefig(buffer, format = 'png', dpi = dpi,
                    pil_kwargs = {'compress_level': 1})
        return buffer.getvalue()

    slider = IntSlider(value = 0, min = 0, max = 180, step = 1, description = 'angle',
                       continuous_update = True,
                       layout = Layout(width = f"{width_in * dpi}px"))

    ## Fixing the displayed size stops the browser from reflowing the layout while the
    ## next frame decodes, which is the other half of the flicker.
    image = Image(value = png_bytes(0), format = 'png',
                  layout = Layout(width = f"{width_in * dpi}px",
                                  height = f"{height_in * dpi}px"))

    def render(change):
        image.value = png_bytes(change['new'])

    slider.observe(render, names = 'value')
    return VBox([slider, image])