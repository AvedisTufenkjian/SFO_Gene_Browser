# app.py

from flask import Flask, render_template, request, jsonify, session
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import umap.umap_ as umap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import umap
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import codecs
#from werkzeug.middleware.profiler import ProfilerMiddleware

app = Flask(__name__)
app.secret_key = ''  

# Load csv files into pandas dataframes
#df1 = pd.read_csv('data/SFO_cont_neurons.csv', index_col=0).T
#df2 = pd.read_csv('data/SFO_cont_cells.csv', index_col=0).T
#df3 = pd.read_csv('data/SFO_hyna_neurons.csv', index_col=0).T  
#df4 = pd.read_csv('data/SFO_hyna_cells.csv', index_col=0).T
#df5 = pd.read_csv('data/SFO_hvol_neurons.csv', index_col=0).T  
#df6 = pd.read_csv('data/SFO_hvol_cells.csv', index_col=0).T  
#df7 = pd.read_csv('data/SFO_wd_neurons.csv', index_col=0).T
#df8 = pd.read_csv('data/SFO_wd_cells.csv', index_col=0).T


##########################################Neurons###############################################################

# Create AnnData objects for each dataset
#adata_1 = sc.AnnData(df1)
#adata_2 = sc.AnnData(df3)
#adata_3 = sc.AnnData(df5)
#adata_4 = sc.AnnData(df7)

# Assign condition labels
#adata_1.obs['condition'] = 'Sated'
#adata_2.obs['condition'] = 'Osmotic'
#adata_3.obs['condition'] = 'Hypovolaemic'
#adata_4.obs['condition'] = 'Water Deprivation'

# Merge datasets Neurons
#adata_merged = adata_1.concatenate(adata_2, adata_3, adata_4, join='outer')
#adata_merged.var_names_make_unique()


# Preprocessing neurons/all?
#sc.pp.normalize_total(adata_merged, target_sum=1e4)
#sc.pp.log1p(adata_merged)
#sc.pp.highly_variable_genes(adata_merged, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.scale(adata_merged, max_value=10)
#sc.tl.pca(adata_merged, svd_solver='arpack')
#sc.pp.neighbors(adata_merged, n_neighbors=10, n_pcs=20)
#sc.tl.umap(adata_merged)
#sc.tl.leiden(adata_merged, key_added='clusters', resolution=2) #ADD IN HIGHER LATER

#adata_merged.obs['clusters'] =adata_merged.obs['clusters'].astype('category')
# create dictionary to map cluster to annotation label
#cluster2annotation = {
#     '0': 'Gaba3-Ptprt',
#     '1': 'Glut4-Il1rapl2',
#     '2': 'Gaba1-Avpr1a', #?
#     '3': 'Glut5-Rxfp3',
#     '4': 'Gaba2-Pnoc',
#     '5': 'Glut2-Calcr',
#     '6': 'Gaba2-Pnoc',
#     '7': 'Glut1-Htr7',
#     '8': 'Glut4-Il1rapl2',
#     '9': 'Gaba1-Avpr1a',
#     '10':'Glut3-Epha5',
#     '11':'Gaba2-Pnoc',
#     '12':'Glut5-Rxfp3',
#     '13':'Glut3-Epha5',
#     '14':'Gaba1-Avpr1a',
#     '15': 'Gaba3-Ptprt',
#     '16': 'Glut4-Il1rapl2',
#     '17': 'Glut4-Il1rapl2',
#    '18':'Gaba2-Pnoc',
#     '19':'Gaba1-Avpr1a',
#     '20':'Gaba1-Avpr1a',
#     '21':'Gaba1-Avpr1a',
#     '22': 'Gaba1-Avpr1a', #?
#     '23':'Gaba1-Avpr1a',
#     '24':'Gaba1-Avpr1a',
#     '25':'Glut2-Calcr',
#     '26':'Glut4-Il1rapl2',
#     '27':'Glut1-Htr7', 
#     '28':'Glut4-Il1rapl2',
#     '29':'Glut2-Calcr',
#     '30':'Gaba2-Pnoc',
#     '31':'Gaba1-Avpr1a',
#     '32':'Gaba2-Pnoc',
#     '33':'Glut5-Rxfp3', #7Glut4-Il1rapl2
#     '34':'Glut4-Il1rapl2',
#     '35':'Glut2-Calcr',
#     '36':'Gaba1-Avpr1a',
#     '37':'37',
#     '38':'38',
#     '39':'39',
#     '40':'40',
#     '41':'41',
#     '42':'42',
#     '43':'43',
#     '44':'44',
#     '45':'45'
#}
# add new .obs column called cell type by mapping clusters to annotation using pandas map function
#adata_merged.obs['Cell Type'] = adata_merged.obs['clusters'].map(cluster2annotation).astype('category')

# Harmony integration
#sc.external.pp.harmony_integrate(adata_merged, key='condition')

# Compute UMAP for the entire dataset
#umap_params = {
#    'n_neighbors': 15,
#    'min_dist': 0.5,
#    'n_components': 2,
#    'metric': 'euclidean'
#}

#umap_embedding = umap.UMAP(**umap_params).fit_transform(adata_merged.obsm['X_pca_harmony'])
#adata_merged.obsm['X_umap'] = umap_embedding

adata_merged = sc.read('data/dataconcatenatedharmonized2_neurons.h5ad')
#adata_merged3d = sc.read('data/dataconcatenatedharmonized2_3d-umap_neurons.h5ad')

# Define manual colors for each cluster
#cluster_colors_neurons = {
#    'Gaba1-Avpr1a': '#1F77B4',
#    'Gaba2-Pnoc': '#FE7F0E',
#   'Gaba3-Ptprt': '#279E68',
#    'Glut1-Htr7':'#D62728',
#    'Glut2-Calcr':'#AA40FB',
#    'Glut3-Epha5':'#8B564B',
#    'Glut4-Il1rapl2':'#E377C2',
#    'Glut5-Rxfp3':'#B5BD61'
#}

##########################################Neurons###############################################################


##########################################Cells#################################################################

# Create AnnData objects for each dataset
#adata_5 = sc.AnnData(df2)
#adata_6 = sc.AnnData(df4)
#adata_7 = sc.AnnData(df6)
#adata_8 = sc.AnnData(df8)

# Assign condition labels
#adata_5.obs['condition'] = 'Sated'
#adata_6.obs['condition'] = 'Osmotic'
#adata_7.obs['condition'] = 'Hypovolaemic'
#adata_8.obs['condition'] = 'Water Deprivation'

# Merge datasets Neurons
#adata_merged2 = adata_5.concatenate(adata_6, adata_7, adata_8, join='outer')
#adata_merged2.var_names_make_unique()

# Preprocessing /all?
#sc.pp.normalize_total(adata_merged2, target_sum=1e4)
#sc.pp.log1p(adata_merged2)
#sc.pp.highly_variable_genes(adata_merged2, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.scale(adata_merged2, max_value=10)
#sc.tl.pca(adata_merged2, svd_solver='arpack')
#sc.pp.neighbors(adata_merged2, n_neighbors=10, n_pcs=20)
#sc.tl.leiden(adata_merged, resolution=0.25)
#sc.tl.umap(adata_merged2)
#sc.tl.leiden(adata_merged2, key_added='clusters', resolution=2) #ADD IN HIGHER LATER

#adata_merged2.obs['clusters'] =adata_merged2.obs['clusters'].astype('category')
#create dictionary to map cluster to annotation label
#cluster2annotation = {
#     '0': 'Inh. Neurons',
#     '1': 'LT Astrocytes',
#     '2': 'Ependymal cell',
#     '3': 'LT Astrocytes',
#     '4': 'LT Astrocytes',
#     '5': 'Astrocytes',
#     '6': 'Inh. Neurons',
#     '7': 'LT Astrocytes',
#     '8': 'LT Astrocytes',
#     '9': 'LT Astrocytes',
#     '10':'Exc. Neurons',
#     '11':'Inh. Neurons',
#     '12':'LT Astrocytes',
#     '13':'LT Astrocytes',
#     '14':'Exc. Neurons',
#     '15':'LT Astrocytes',
#     '16':'LT Astrocytes',
#     '17':'Exc. Neurons',
#     '18':'Inh. Neurons',
#     '19':'LT Astrocytes',
#     '20':'LT Astrocytes',
#     '21':'Exc. Neurons',
#     '22':'LT Astrocytes',
#     '23':'Exc. Neurons',
#     '24':'LT Endo Cells',
#     '25':'LT Astrocytes',
#     '26':'Fibroblasts',
#     '27':'Microglia',
#     '28':'LT Endo Cells',
#     '29':'Inh. Neurons',
#     '30':'Endo Cells',
#     '31':'LT Astrocytes',
#     '32':'Inh. Neurons',
#     '33':'Exc. Neurons',
#     '34':'LT Astrocytes',
#     '35':'LT Astrocytes',
#     '36':'Pericytes', #?what looks better
#     '37':'Microglia',
#     '38':'Exc. Neurons',
#     '39':'Fibroblasts',
#     '40':'Pericytes',
#     '41':'Microglia',
#     '42':'Ependymal cell',
#     '43':'Oligodendrocytes',
#     '44':'Ependymal cell',
#     '45':'VSMCs' #pericytes?
#}

# add new .obs column called cell type by mapping clusters to annotation using pandas map function
#adata_merged2.obs['Cell Type'] = adata_merged2.obs['clusters'].map(cluster2annotation).astype('category')

# Harmony integration
#sc.external.pp.harmony_integrate(adata_merged2, key='condition')

# Compute UMAP for the entire dataset
#umap_params = {
#    'n_neighbors': 15,
#    'min_dist': 0.5,
#    'n_components': 2,
#    'metric': 'euclidean'
#}

#umap_embedding = umap.UMAP(**umap_params).fit_transform(adata_merged2.obsm['X_pca_harmony'])
#adata_merged2.obsm['X_umap'] = umap_embedding

adata_merged2 = sc.read('data/dataconcatenatedharmonized2_cells.h5ad')
#adata_merged23d = sc.read('data/dataconcatenatedharmonized2_3d-umap_cells.h5ad')


# Define manual colors for each cluster
#cluster_colors_cells = {
#   'Astrocytes': '#1F77B4',
#    'Endo Cells': '#FE7F0E',
#    'Ependymal cell': '#279E68',
#    'Exc. Neurons':'#D62728',
#    'Fibroblasts':'#AA40FB',
#    'Inh. Neurons':'#8B564B',
#    'LT Astrocytes':'#E377C2',
#    'LT Endo Cells':'#B5BD61',
#    'Microglia':'#17BECF',
#    'Oligodendrocytes':'#AEC7E8',
#    'Pericytes':'#FEBB78',
#    'VSMCs':'#98DF8A'
#}

##########################################Cells#################################################################


def generate_scanpy_plots(gene_name, adata, selected_file):
    # Create a copy of the adata object for the 3D UMAP
    # adata_3d = adata.copy()

    # Compute 2D UMAP for the original adata object- not necessary ?
    # sc.tl.umap(adata, n_components=2)

    # Compute 3D UMAP for the copied adata object
    # sc.tl.umap(adata_3d, n_components=3)

    expression_values = adata.obs_vector(gene_name)
    vmin = np.min(expression_values)
    vmax = np.max(expression_values)

    fig = plt.figure(figsize=(20.3, 13))  # 23,17

    # Adjust height_ratios to include an extra row for white space
    gs = gridspec.GridSpec(3, 1, height_ratios=[.95, 0.05, 2.25])  # Added 0.2 for white space 0.75, 0.05, 2.25

    # Create a smaller GridSpecFromSubplotSpec for the first row
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[0])

    conditions = ['Sated', 'Osmotic', 'Hypovolaemic', 'Water Deprivation', 'Sodium Depletion']

    for i, condition in enumerate(conditions):
        ax = fig.add_subplot(gs1[0, i])
        adata_condition = adata[adata.obs['condition'] == condition]
        sc.pl.umap(adata_condition, color=gene_name, legend_loc='on data', vmin=vmin, vmax=vmax,
                   title=f'{gene_name}\n under {condition} condition', ax=ax, show=False)

    # Create a smaller GridSpecFromSubplotSpec for the second and third rows (ignoring the middle row for white space)
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 5, subplot_spec=gs[2], hspace=0.6)  # Keep hspace as before

    # Determine global min and max across all conditions
    global_min, global_max = np.inf, -np.inf
    for condition in conditions:
        adata_condition = adata[adata.obs['condition'] == condition]
        expression_values = adata_condition.obs_vector(gene_name)
        vmin, vmax = np.min(expression_values), np.max(expression_values)
        global_min, global_max = min(global_min, vmin), max(global_max, vmax)

    for i, condition in enumerate(conditions):
        ax = fig.add_subplot(gs2[0, i])
        adata_condition = adata[adata.obs['condition'] == condition]
        sc.pl.violin(adata_condition, keys=gene_name, rotation=90, groupby='Cell Type', ax=ax, show=False)
        ax.set_ylim([global_min, global_max])  # Set y-axis limits here

    ax = fig.add_subplot(gs2[1, 1])
    sc.pl.dotplot(adata, var_names=gene_name, groupby='condition', dendrogram=True, ax=ax, show=False)

    if adata is adata_merged:
        title = 'Neuron Types in the SFO'
    elif adata is adata_merged2:
        title = 'Cell Types in the SFO'

    ax = fig.add_subplot(gs2[1, 2])
    sc.pl.umap(adata, color='Cell Type', frameon=False, title=title, ax=ax, show=False)

    # Use tight_layout to automatically adjust subplot parameters
    plt.tight_layout()

    # Before saving the figure, adjust the subplots to remove white space
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)  # Adjust the margins to reduce white space

    # Save the Matplotlib figure to a BytesIO object
    bio_matplotlib = BytesIO()
    plt.savefig(bio_matplotlib, format='png', bbox_inches='tight')
    bio_matplotlib.seek(0)

    # Convert the Matplotlib figure to a base64 string
    img_base64 = base64.b64encode(bio_matplotlib.getvalue()).decode('utf-8')

    plt.close()



    # Determine which cluster colors to use based on the type of adata 
#    if adata is adata_merged: 
#        cluster_colors = cluster_colors_neurons 
#    elif adata is adata_merged2: 
#        cluster_colors = cluster_colors_cells 
    #else: # Default to some colors jata object doesn't match known types 

#    fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

#    for cluster, color in cluster_colors.items():
#        cluster_data = adata_3d[adata_3d.obs['Cell Type'] == cluster]  

        # Define the size of the dots based on a condition related to the adata object
#        if selected_file == 'Neurons':
#            size = 2
#        elif selected_file == 'Cells':
#            size = 1

        # Title
#        if selected_file == 'Neurons':
#            title2 = '3D UMAP of Neuron Types in the SFO<br><sup>drag to pan, scroll to zoom, double click on legend to isolate single neuron class in the SFO</sup>'
#        elif selected_file == 'Cells':
#            title2 = '3D UMAP of Cell Types in the SFO<br><sup>drag to pan, scroll to zoom, double click on legend to isolate single cell class in the SFO</sup>'


#        fig.add_trace(go.Scatter3d(
#            x=cluster_data.obsm['X_umap'][:, 0],
#            y=cluster_data.obsm['X_umap'][:, 1],
#            z=cluster_data.obsm['X_umap'][:, 2],
#            mode='markers',
#            marker=dict(
#                size=size,
#                color=color,
#                opacity=0.8,
#            ),
#            name=cluster
#        ))

#    fig.update_layout(
#        title=title2,
#        scene=dict(
#            xaxis=dict(title='UMAP 1'),
#            yaxis=dict(title='UMAP 2'),
#            zaxis=dict(title='UMAP 3')
#        ),
#        autosize=False,
#        width=800,  # adjust as needed
#        height=640,  # adjust as needed
#    )

#little if statement
    if selected_file == 'Neurons':
            html_file_path = "data/3d-umap_neurons_plotly_plot2.html"
    elif selected_file == 'Cells':
            html_file_path = "data/3d-umap_cells_plotly_plot2.html"
    
    # Convert Plotly figure to HTML div string
    #plot_html = pio.to_html(fig, full_html=False, config={'displaylogo': False})
#    html_file_path = "data/3d-umap_neurons_plotly_plot2.html"
    with codecs.open(html_file_path, 'r', 'utf-8') as f:
        plot_html = f.read()
    #plot_html = 'data/3d-umap_neurons_plotly_plot.html'

    return img_base64, plot_html



# Function to get autofill suggestions based on user input and adataa object
def get_autofill_suggestions(prefix, adata):
    # adata.var_names is list of gene names
    suggestions = [gene for gene in adata.var_names if gene.lower().startswith(prefix)]

    # Limit number of suggestions on search autofill results. if needed
    return suggestions[:10]

# Route for selecting csv file
@app.route('/select_file', methods=['POST'])
def select_file():
    selected_file = request.form['csv_file']
    session['selected_file'] = selected_file

    # Load csv file based on selected option
    if selected_file == 'Neurons':
        adata = adata_merged
    elif selected_file == 'Cells':
        adata = adata_merged2 #CHANGE ONCE CELL DATASET IS ACTIVE and the rest below
    else:
        # Default to adtat1 if nothing is selected
        adata = adata_merged

    return render_template('index.html.html', selected_file=selected_file)

# Route for main page
@app.route('/')
def index():
    return render_template('file_selection.html.html')

# Route for plotting single genes
@app.route('/search', methods=['POST'])
def search():
    search_value = request.form['search']
    selected_file = session.get('selected_file')

    # Determine which adata object to use based on selected file
    if selected_file == 'Neurons':
        adata = adata_merged
    elif selected_file == 'Cells':
        adata = adata_merged2
    else:
        # Default to adata1 object else
        adata = adata_merged2

    # Generate scanpy plots for the specified gene and adata object
    img_base64, plot_html = generate_scanpy_plots(search_value, adata, selected_file) #working# img_base64, plot_html = generate_scanpy_plots(search_value, adata, selected_file)


    return render_template('index.html.html', selected_file=selected_file, search_value=search_value, img_base64=img_base64, plot_html=plot_html) #working# return render_template('index.html.html', selected_file=selected_file, search_value=search_value, img_base64=img_base64, plot_html=plot_html)




# Route for plotting multiple genes
#@app.route('/plot_multiple_genes', methods=['GET', 'POST'])
#def plot_multiple_genes():
    selected_file = session.get('selected_file')

    # Determine adata object to use based on selected library
    if selected_file == 'Neurons':
        adata = adata_merged
    elif selected_file == 'Cells':
        adata = adata_merged2
    else:
        # Default to adata1 if nothing is selected
        adata = adata_merged
    if request.method == 'POST':
        # Get list of genes from user input
        gene_list = request.form['gene_list'].split(',')

        # Generate scanpy plots for multiple genes
        img_base64 = generate_multiple_genes_plot(gene_list, adata)

        return render_template('plot_multiple_genes.html.html', selected_file=selected_file, gene_list=gene_list, img_base64=img_base64)

    return render_template('plot_multiple_genes.html.html', selected_file=selected_file)

# Function to generate scanpy plots for multiple genes
#def generate_multiple_genes_plot(gene_names, adata):
    # Create subplots - original fig size (14,5), increased to avoid overlap, tight.plt not working??
    fig, axs = plt.subplots(1, 3, figsize=(25, 5))

    # Dot plot for multiple genes
    sc.pl.dotplot(adata, var_names=gene_names, groupby='Cell Type', dendrogram=True, ax=axs[2], show=False)

    # Stacked violin plot for multiple genes
    sc.pl.stacked_violin(adata, var_names=gene_names, groupby='Cell Type', swap_axes=False, ax=axs[1], show=False)

    # Matrix plot for multiple genes
    sc.pl.matrixplot(adata, var_names=gene_names, groupby='Cell Type', swap_axes=False, ax=axs[0], show=False)

    # Save the plots to a BytesIO object
    bio = BytesIO()
    plt.savefig(bio, format='png')
    bio.seek(0)

    # Encode the image to base64 for embedding in HTML
    img_base64 = base64.b64encode(bio.getvalue()).decode('utf-8')

    # Close plot
    plt.close()

    return img_base64

# Route to get autofill suggestions
@app.route('/autofill', methods=['POST'])
def autofill():
    search_prefix = request.form['search_value'].lower()
    selected_file = session.get('selected_file')

    # Determine which adata object to use based off selected file -  necessary
    if selected_file == 'Neurons':
        adata = adata_merged
    elif selected_file == 'Cells':
        adata = adata_merged
    else:
        # Default to the adata1 if no file is selected
        adata = adata

    suggestions = get_autofill_suggestions(search_prefix, adata)
    return jsonify(suggestions)

#app.wsgi_app = ProfilerMiddleware(app.wsgi_app)


@app.route('/Jerboa Brain Atlas')
def tool2():
    return render_template('tool2.html')



if __name__ == '__main__':
    app.run(port=80) #Port 80 before export !!
