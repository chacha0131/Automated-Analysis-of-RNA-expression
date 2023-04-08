import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import ConvexHull

files = next(os.walk('/Users/emilycha/Desktop/imgs_bagot/PVT/out/quant/'))[1]

for file in files:
    print(file)
    path = os.path.join('/Users/emilycha/Desktop/imgs_bagot/PVT/out/quant/', file)
    mlcls = pd.read_csv(os.path.join(path, file + '.csv'))
    cells = pd.read_csv(os.path.join(path, file + '_cell_stats.csv'))
    n_cells = cells.shape[0]
    mlcls_noise = mlcls[mlcls['cell'] == 0]
    mlcls_cells = mlcls[mlcls['cell'] != 0]
    mlcls_per_cell = mlcls_cells.groupby(['gene','cell']).size().reset_index(name='counts')
    mlcls_cols = {'vGlut_denoised':'paleturquoise','vGAT_denoised':'palegreen','Sdk1_denoised':'pink'}
    cells_cols = {'vgat+/vglut-/sdk1-':'palegreen','vgat+/vglut-/sdk1+':'green','vgat+/vglut+/sdk1-':'cornflowerblue',
                    'vgat-/vglut+/sdk1-':'paleturquoise','vgat-/vglut+/sdk1+':'turquoise','vgat+/vglut+/sdk1+':'cornflowerblue',
                    'vgat-/vglut-/sdk1+':'orangered', 'none':'gold'}
    cells['type'] = 'cell_type'
    cells['vgat'] = 0
    cells['vglut'] = 0
    cells['sdk1'] = 0
    for cell in cells.cell:
        try:
            vgat = mlcls_per_cell.loc[(mlcls_per_cell['cell'] == cell) & (mlcls_per_cell['gene'] == 'vGAT_denoised')].counts.reset_index(drop=True)[0]
        except KeyError:
            vgat = 0
        try: 
            vglut = mlcls_per_cell.loc[(mlcls_per_cell['cell'] == cell) & (mlcls_per_cell['gene'] == 'vGlut_denoised')].counts.reset_index(drop=True)[0]
        except KeyError:
            vglut = 0
        try:
            sdk1 = mlcls_per_cell.loc[(mlcls_per_cell['cell'] == cell) & (mlcls_per_cell['gene'] == 'Sdk1_denoised')].counts.reset_index(drop=True)[0]
        except KeyError:
            sdk1 = 0
        if sdk1 < 2:
            if (vgat >= 10) & (vglut < 13):
                cell_type = 'vgat+/vglut-/sdk1-'
            elif (vgat < 10) & (vglut >= 13):
                cell_type = 'vgat-/vglut+/sdk1-'
            elif (vgat >= 10) & (vglut >= 13):
                if vgat/(vglut+vgat) > (2/3):
                    cell_type = 'vgat+/vglut-/sdk1-'
                elif vglut/(vglut+vgat) > (2/3):
                    cell_type = 'vgat-/vglut+/sdk1-'
                else:
                    cell_type = 'vgat+/vglut+/sdk1-'
            elif (vgat < 10) & (vglut < 13):
                cell_type = 'none'
        elif sdk1 >= 2:
            if (vgat >= 10) & (vglut < 13):
                cell_type = 'vgat+/vglut-/sdk1+'
            elif (vgat < 10) & (vglut >= 13):
                cell_type = 'vgat-/vglut+/sdk1+'
            elif (vgat >= 10) & (vglut >= 13):
                if vgat/(vglut+vgat) > (2/3):
                    cell_type = 'vgat+/vglut-/sdk1+'
                elif vglut/(vglut+vgat) > (2/3):
                    cell_type = 'vgat-/vglut+/sdk1+'
                else:
                    cell_type = 'vgat+/vglut+/sdk1+'
            elif (vgat < 10) & (vglut < 13):
                if sdk1 >= 3:
                    cell_type = 'vgat-/vglut-/sdk1+'
                else:
                    cell_type = 'none'
        cells.loc[cells.cell == cell,'type'] = cell_type
        cells.loc[cells.cell == cell,'vgat'] = vgat
        cells.loc[cells.cell == cell,'vglut'] = vglut
        cells.loc[cells.cell == cell,'sdk1'] = sdk1
    cells.to_csv(path + '/cells_quant.csv')
    
    
    cell_types = pd.Series(cells.type)
    cell_types_count = cell_types.value_counts() 
    pie_chart = cell_types_count.plot(kind='pie')
    pie_chart.figure.savefig(path + '/type_chart.png')
    plt.close(pie_chart.get_figure())
    im = plt.imread('/Users/emilycha/Desktop/imgs_bagot/PVT/imgs/preproc/nucl_raw/' + file + '_nuclei.tif')
    fig = plt.figure()
    px = 1/plt.rcParams['figure.dpi']
    fig.set_size_inches(im.shape[1]*px,im.shape[0]*px)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    implot = plt.imshow(im, cmap = 'gray')
    plt.scatter(mlcls_cells.x, mlcls_cells.y, c = mlcls_cells['gene'].map(mlcls_cols), s = 2)
    for cell in cells.cell:
        if cell % 200 == 0:
            print(str(cell) + ' / ' + str(n_cells))
        cell_mlcls = mlcls_cells.loc[mlcls_cells['cell'] == cell].reset_index()
        cell_points = np.array([cell_mlcls.x,cell_mlcls.y]).transpose().astype(int)
        if cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0] in ['none']:
            continue
        try:
            cell_hull = ConvexHull(cell_points)
            cell_col = cells_cols[cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0]]
            for simplex in cell_hull.simplices:
                plt.plot(cell_points[simplex, 0], cell_points[simplex, 1], c = cell_col)
        except Exception as e:
            continue
    plt.savefig(path + '/all_cells.png')
    plt.close()
    im = plt.imread('/Users/emilycha/Desktop/imgs_bagot/PVT/imgs/preproc/nucl_raw/' + file + '_nuclei.tif')
    fig = plt.figure()
    px = 1/plt.rcParams['figure.dpi']
    fig.set_size_inches(im.shape[1]*px,im.shape[0]*px)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    implot = plt.imshow(im, cmap = 'gray')
    plt.scatter(mlcls_cells.x, mlcls_cells.y, c = mlcls_cells['gene'].map(mlcls_cols), s = 2)
    for cell in cells.cell:
        if cell % 200 == 0:
            print(str(cell) + ' / ' + str(n_cells))
        cell_mlcls = mlcls_cells.loc[mlcls_cells['cell'] == cell].reset_index()
        cell_points = np.array([cell_mlcls.x,cell_mlcls.y]).transpose().astype(int)
        if cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0] in ['none','vgat+/vglut-/sdk1-','vgat-/vglut+/sdk1-','vgat+/vglut+/sdk1-']:
            continue
        try:
            cell_hull = ConvexHull(cell_points)
            cell_col = cells_cols[cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0]]
            for simplex in cell_hull.simplices:
                plt.plot(cell_points[simplex, 0], cell_points[simplex, 1], c = cell_col)
        except Exception as e:
            continue
    plt.savefig(path + '/pos_cells.png')
    plt.close()
    im = plt.imread('/Users/emilycha/Desktop/imgs_bagot/PVT/imgs/preproc/nucl_raw/' + file + '_nuclei.tif')
    fig = plt.figure()
    px = 1/plt.rcParams['figure.dpi']
    fig.set_size_inches(im.shape[1]*px,im.shape[0]*px)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    implot = plt.imshow(im, cmap = 'gray')
    plt.scatter(mlcls_cells.x, mlcls_cells.y, c = mlcls_cells['gene'].map(mlcls_cols), s = 2)
    for cell in cells.cell:
        if cell % 200 == 0:
            print(str(cell) + ' / ' + str(n_cells))
        cell_mlcls = mlcls_cells.loc[mlcls_cells['cell'] == cell].reset_index()
        cell_points = np.array([cell_mlcls.x,cell_mlcls.y]).transpose().astype(int)
        if cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0] not in ['none']:
            continue
        try:
            cell_hull = ConvexHull(cell_points)
            cell_col = cells_cols[cells.loc[cells['cell'] == cell].type.reset_index(drop=True)[0]]
            for simplex in cell_hull.simplices:
                plt.plot(cell_points[simplex, 0], cell_points[simplex, 1], c = cell_col)
        except Exception as e:
            continue
    plt.savefig(path + '/noise_cells.png')
    plt.close()
    

