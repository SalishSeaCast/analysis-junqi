# In this script, we will generate the plots of sectors used in the appendix of our report.
# Why are we using script instead of a notebllok?
# It's because I'm in Bamfield and I cannot connect to our Jupyter core through ssh.

import netCDF4 as nc
import numpy as np
import datetime
import os
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm
import string

# --- 1. 配置与路径设置 ---
BASE_DIR = '/results2/SalishSea/nowcast-green.201905/'
FNAME_HEAD = 'SalishSea_1d_' 
FNAME_TAIL = '_ptrc_T.nc'
NEW_INDIR_RESULTS = '/home/jqiu/Programing/Projects/analysis-junqi/Diatom_vs_Flagellate_Report/Results_Section_Biology/Results_Section_pdf/'

if not os.path.exists(NEW_INDIR_RESULTS):
    os.makedirs(NEW_INDIR_RESULTS)

TARGETS_DEFINITIONS = {
    'Point_1': {'value': (49.0, -123.25), 'label': 'Section#1'},
    'Point_2': {'value': (49.3, -124.0),  'label': 'Section#2'},
    'Point_3': {'value': (49.9, -124.8),  'label': 'Section#3'},
}

SECTION_CONFIG = {
    'Point_1': {'left_offset': 35, 'right_offset': 22},
    'Point_2': {'left_offset': 30, 'right_offset': 70},
    'Point_3': {'left_offset': 45, 'right_offset': 25},
}

year_target = [2011, 2018]
month_target = 7
DAY_PRE = 9   
DAY_POST = 17 
DEPTH_LIMIT = 50

# --- 2. 辅助函数 ---
def find_nearest_indices_2d(lon_2d, lat_2d, target_lon, target_lat):
    dist_sq = (lon_2d - target_lon)**2 + (lat_2d - target_lat)**2
    return np.unravel_index(dist_sq.argmin(), lon_2d.shape)

def load_bio_data_refined(date_obj, y_idx, x_idx, config, lon_2d):
    file_path = os.path.join(BASE_DIR, date_obj.strftime('%d%b%y').lower(), 
                             f"{FNAME_HEAD}{date_obj.strftime('%Y%m%d')}_{date_obj.strftime('%Y%m%d')}{FNAME_TAIL}")
    if not os.path.exists(file_path): return None
    
    with nc.Dataset(file_path, 'r') as ncfile:
        z_var = ncfile.variables['deptht'][:]
        z_idx_limit = np.abs(z_var - DEPTH_LIMIT).argmin()
        x_start = max(0, x_idx - config['left_offset'])
        x_end = min(lon_2d.shape[1], x_idx + config['right_offset'])
        x_slice = slice(x_start, x_end)
        
        diat = ncfile.variables['diatoms'][0, :z_idx_limit+1, y_idx, x_slice]
        flag = ncfile.variables['flagellates'][0, :z_idx_limit+1, y_idx, x_slice]
        lon_sec = lon_2d[y_idx, x_slice]
        return diat, flag, lon_sec, z_var[:z_idx_limit+1]

# --- 3. 核心绘图程序 ---

def main():
    ref_path = os.path.join(BASE_DIR, '08jul11/SalishSea_1d_20110708_20110708_ptrc_T.nc')
    with nc.Dataset(ref_path, 'r') as nf:
        lon_2d = nf.variables['nav_lon'][:]
        lat_2d = nf.variables['nav_lat'][:]

    print("Loading data...")
    data_cache = {y: {pid: {} for pid in TARGETS_DEFINITIONS} for y in year_target}
    
    for year in year_target:
        for p_id, p_def in TARGETS_DEFINITIONS.items():
            y_idx, x_idx = find_nearest_indices_2d(lon_2d, lat_2d, p_def['value'][1], p_def['value'][0])
            res_pre = load_bio_data_refined(datetime.date(year, month_target, DAY_PRE), y_idx, x_idx, SECTION_CONFIG[p_id], lon_2d)
            res_post = load_bio_data_refined(datetime.date(year, month_target, DAY_POST), y_idx, x_idx, SECTION_CONFIG[p_id], lon_2d)
            
            if res_pre and res_post:
                data_cache[year][p_id]['lon'] = res_pre[2]
                data_cache[year][p_id]['depth'] = res_pre[3]
                data_cache[year][p_id]['diat_pre'], data_cache[year][p_id]['flag_pre'] = res_pre[0], res_pre[1]
                data_cache[year][p_id]['diat_post'], data_cache[year][p_id]['flag_post'] = res_post[0], res_post[1]

    # ------------------------------------------------------------------
    # 图 1：差值图 (July 17 - July 9)  [4行 x 3列] (保持不变)
    # ------------------------------------------------------------------
    print("Plotting Difference Figure...")
    fig1, axes1 = plt.subplots(4, 3, figsize=(16, 14), sharex='col', sharey='row', gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    fig1.patch.set_facecolor('white')
    
    row_configs_diff = [
        {'year': 2011, 'var': 'diat', 'title': '2011 Diatom Diff', 'vmax': 10},
        {'year': 2018, 'var': 'diat', 'title': '2018 Diatom Diff', 'vmax': 10},
        {'year': 2011, 'var': 'flag', 'title': '2011 Flagellate Diff', 'vmax': 5},
        {'year': 2018, 'var': 'flag', 'title': '2018 Flagellate Diff', 'vmax': 5}
    ]
    p_keys = list(TARGETS_DEFINITIONS.keys())
    abc_labels = list(string.ascii_lowercase)
    
    for r, r_cfg in enumerate(row_configs_diff):
        y = r_cfg['year']
        v_type = r_cfg['var']
        for c, p_id in enumerate(p_keys):
            ax = axes1[r, c]
            ax.set_facecolor('darkgray')
            d_post = data_cache[y][p_id][f'{v_type}_post']
            d_pre = data_cache[y][p_id][f'{v_type}_pre']
            diff = d_post - d_pre
            diff_masked = np.ma.masked_where((d_pre <= 0) & (d_post <= 0), diff)
            
            X, Z = np.meshgrid(data_cache[y][p_id]['lon'], data_cache[y][p_id]['depth'])
            cf = ax.pcolormesh(X, Z, diff_masked, cmap='RdBu_r', vmin=-r_cfg['vmax'], vmax=r_cfg['vmax'], shading='auto')
            
            if r == 0: ax.set_title(TARGETS_DEFINITIONS[p_id]['label'], fontsize=14, fontweight='bold')
            if c == 0: ax.set_ylabel(f"{r_cfg['title']}\nDepth (m)", fontsize=12)
            if r == 3: ax.set_xlabel('Longitude ($^{\circ}$E)', fontsize=12)
            
            t = ax.text(0.03, 0.95, f"({abc_labels[r*3 + c]})", transform=ax.transAxes, ha='left', va='top', fontsize=14, fontweight='bold', color='black')
            t.set_path_effects([pe.withStroke(linewidth=2, foreground='white')])
            
            if c == 2 and r in [0, 2]:
                cbar_ax = fig1.add_axes([0.92, 0.53 if r==0 else 0.13, 0.015, 0.35])
                cbar = fig1.colorbar(cf, cax=cbar_ax, extend='both')
                v_name = "Diatom" if r==0 else "Flagellate"
                cbar.set_label(f'$\Delta$ {v_name} ($\mu$mol N L$^{{-1}}$)', fontsize=12)

    axes1[0,0].invert_yaxis()
    fig1.savefig(os.path.join(NEW_INDIR_RESULTS, "Difference_Panel_July17_minus_9.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig1)

    # ------------------------------------------------------------------
    # 图 2：比值图 (Diatom / Flagellate) - 包含风暴前后 [4行 x 3列]
    # ------------------------------------------------------------------
    print("Plotting Ratio Figure (Pre & Post Storm)...")
    fig2, axes2 = plt.subplots(4, 3, figsize=(16, 14), sharex='col', sharey='row', gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    fig2.patch.set_facecolor('white')
    
    # 按照 年份 -> 风暴前 -> 风暴后 的顺序排列
    row_configs_ratio = [
        {'year': 2011, 'state': 'pre',  'title': '2011 Pre (Jul 9)'},
        {'year': 2011, 'state': 'post', 'title': '2011 Post (Jul 17)'},
        {'year': 2018, 'state': 'pre',  'title': '2018 Pre (Jul 9)'},
        {'year': 2018, 'state': 'post', 'title': '2018 Post (Jul 17)'}
    ]
    
    for r, r_cfg in enumerate(row_configs_ratio):
        y = r_cfg['year']
        state = r_cfg['state']
        for c, p_id in enumerate(p_keys):
            ax = axes2[r, c]
            ax.set_facecolor('darkgray')
            
            diat = data_cache[y][p_id][f'diat_{state}']
            flag = data_cache[y][p_id][f'flag_{state}']
            
            ratio = diat / (flag + 1e-6)
            ratio_masked = np.ma.masked_where((diat <= 0.01) | (flag <= 0.01), ratio)
            
            X, Z = np.meshgrid(data_cache[y][p_id]['lon'], data_cache[y][p_id]['depth'])
            
            cf2 = ax.pcolormesh(X, Z, ratio_masked, cmap='PiYG', 
                                norm=LogNorm(vmin=0.1, vmax=10), shading='auto')
            
            if r == 0: ax.set_title(TARGETS_DEFINITIONS[p_id]['label'], fontsize=14, fontweight='bold')
            if c == 0: ax.set_ylabel(f"{r_cfg['title']}\nDepth (m)", fontsize=12)
            if r == 3: ax.set_xlabel('Longitude ($^{\circ}$E)', fontsize=12)
            
            t2 = ax.text(0.03, 0.95, f"({abc_labels[r*3 + c]})", transform=ax.transAxes, 
                         ha='left', va='top', fontsize=14, fontweight='bold', color='black')
            t2.set_path_effects([pe.withStroke(linewidth=2, foreground='white')])

    axes2[0,0].invert_yaxis()
    
    # 统一使用一个大 Colorbar
    cbar_ax2 = fig2.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar2 = fig2.colorbar(cf2, cax=cbar_ax2, extend='both')
    cbar2.set_label('Ratio: Diatom / Flagellate', fontsize=13)
    
    fig2.savefig(os.path.join(NEW_INDIR_RESULTS, "Ratio_Panel_Pre_vs_Post.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig2)

    print(f"All figures successfully generated and saved as PDF in {NEW_INDIR_RESULTS}")

if __name__ == '__main__':
    main()