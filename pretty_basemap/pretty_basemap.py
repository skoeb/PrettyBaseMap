#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Sam Koebrich so that I will never have a lazy map again...
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import os
import shapely

PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

class PrettyStateBaseMap():
    def __init__(self, selected_usstates, dpi=250):
        self.selected_usstates = selected_usstates
        self.usstateshapefiles = os.path.join(PATH, 'geometry', 'cb_2018_us_state_20m', 'cb_2018_us_state_20m.shp')
        self.new_crs = ccrs.Mercator()
        self.dpi = dpi
        
    def _load_shapes(self):
        self.usstates = gpd.read_file(self.usstateshapefiles)
        self.usstates = self.usstates.to_crs(self.new_crs.proj4_init)

    def _filter_islands(self, gdf, percentile=0.05):
        gdf["area"] = gdf['geometry'].area/ 10**6
        p = gdf["area"].quantile(percentile)
        gdf = gdf.loc[gdf['area'] > p]
        return gdf
    
    def _filter_centroid(self, gdf, min_lat=-70):
        gdf['centroid'] = gdf['geometry'].centroid
        gdf['centroid'] = [i.coords[0] for i in gdf['centroid']]
        gdf['centroid_lon'] = [i[0] for i in gdf['centroid']]
        gdf['centroid_lat'] = [i[1] for i in gdf['centroid']]
        
        # --- Filter antarctica ---
        gdf = gdf.loc[gdf['centroid_lat'] > min_lat]
        return gdf
    
    def _apply_filters(self, gdf):
        gdf = self._filter_islands(gdf)
        gdf = self._filter_centroid(gdf)
        return gdf

    def _change_crs(self):
        self.usstates = self.usstates.to_crs(self.new_crs.proj4_init)
        return self
    
    def _filter_shapes(self):
        self.usstates = self.usstates.loc[self.usstates['STUSPS'].isin(self.selected_usstates)]
        
    def _get_bounds(self):
        disolved = self.usstates.unary_union
        bounds = disolved.bounds
        
        # --- Clean bounds ---
        bounds = [round(i, 1) for i in bounds]
        self.bounds = [bounds[0]*1.02, bounds[2]*0.98, bounds[1]
                       * 0.98, bounds[3]*1.02]  # add a bit of buffer

        
    def fig_ax(self, cartopy_features=True):
        self._load_shapes()
        self._filter_shapes()
        self._change_crs()
        self._get_bounds()

        # --- Initialize Figure ---
        fig = plt.figure(dpi=self.dpi)
        ax = fig.add_axes([0, 0, 1, 1], projection=self.new_crs)
        ax.set_extent(self.bounds, self.new_crs)
        

        # --- Make Pretty Background ---
        if cartopy_features:
            stamen_terrain = cimgt.Stamen('terrain-background')
            ax.add_image(stamen_terrain, 5)
            ax.add_feature(cfeature.LAKES)
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.BORDERS)
        
        self.usstates.plot(edgecolor='k', facecolor='None',
                            alpha=1, linewidth=0.5, zorder=10,
                            ax=ax)
        
        # --- Add Source ---
        # fig.text(0, 0.17, 'Figure by NREL', size=7)
        
        ax.axis('off')
        fig.tight_layout()
        
        return fig, ax


class PrettyCountryBaseMap():
    def __init__(self, selected_iso3, selected_usstates=None, dpi=250):
        self.selected_iso3 = selected_iso3
        self.selected_usstates = selected_usstates
        self.countryshapefiles = os.path.join(PATH, 'geometry','ne_110m_admin_0_countries','ne_110m_admin_0_countries.shp')
        self.usstateshapefiles = os.path.join(PATH, 'geometry', 'cb_2018_us_state_20m', 'cb_2018_us_state_20m.shp')
        self.new_crs = ccrs.EuroPP()
        self.old_crs = 'WGS84' #ccrs.PlateCarree()
        self.dpi = dpi
        
    def _load_shapes(self):
        self.country = gpd.read_file(self.countryshapefiles)
        self.country = self.country.to_crs(self.new_crs.proj4_init)

        if self.selected_iso3 == 'USA':
            self.usstates = gpd.read_file(self.usstateshapefiles)
            self.usstates = self.usstates.to_crs(self.new_crs.proj4_init)

    def _filter_islands(self, gdf, percentile=0.05):
        gdf["area"] = gdf['geometry'].area/ 10**6
        p = gdf["area"].quantile(percentile)
        gdf = gdf.loc[gdf['area'] > p]
        return gdf
    
    def _filter_centroid(self, gdf, min_lat=-70):
        gdf['centroid'] = gdf['geometry'].centroid
        gdf['centroid'] = [i.coords[0] for i in gdf['centroid']]
        gdf['centroid_lon'] = [i[0] for i in gdf['centroid']]
        gdf['centroid_lat'] = [i[1] for i in gdf['centroid']]
        
        # --- Filter antarctica ---
        gdf = gdf.loc[gdf['centroid_lat'] > min_lat]
        return gdf
    
    def _apply_filters(self, gdf):
        gdf = self._filter_islands(gdf)
        gdf = self._filter_centroid(gdf)
        return gdf

    def _change_crs(self):
        self.country = self.country.to_crs(self.new_crs.proj4_init)

        if self.selected_iso3 == 'USA':
            self.usstates = self.usstates.to_crs(self.new_crs.proj4_init)
    
    def _filter_shapes(self):
        self.country = self._apply_filters(self.country)
        self.country = self.country.loc[self.country['ISO_A3'] == self.selected_iso3]

        if (self.selected_usstates !=None) & (self.selected_iso3 == 'USA'):
            self.usstates = self.usstates.loc[self.usstates['STUSPS'].isin(self.selected_usstates)]
        
    def _get_bounds(self):
        assert len(self.country) == 1

        # --- Get the geometry of the country ---
        country_geo = self.country['geometry'].item()
        
        # --- Only take the biggest polygon (i.e. no alaska or hawaii for US) ---
        if isinstance(country_geo, shapely.geometry.multipolygon.MultiPolygon):
            areas = [p.area for p in country_geo]
            biggest_area = max(areas)
            country_geo = [p for p in country_geo if p.area == biggest_area][0]
        
        # --- Get bounds of largest polygon ---
        bounds = country_geo.bounds
        
        # --- Clean bounds ---
        bounds = [round(i, 1) for i in bounds]
        self.bounds = [bounds[0]*1, bounds[2]*1, bounds[1]*1, bounds[3]*1]
        
    def fig_ax(self, cartopy_features=True, us_states=True):
        self._load_shapes()
        self._filter_shapes()
        self._change_crs()
        self._get_bounds()

        # --- Initialize Figure ---
        fig = plt.figure(dpi=self.dpi)
        ax = fig.add_axes([0, 0, 1, 1], projection=self.new_crs)
        ax.set_extent(self.bounds, self.new_crs)

        # --- Make Pretty Background ---
        if cartopy_features:
            stamen_terrain = cimgt.Stamen('terrain-background')
            ax.add_image(stamen_terrain, 8)
            # ax.add_feature(cfeature.LAKES)
            # ax.add_feature(cfeature.OCEAN)
            # ax.add_feature(cfeature.BORDERS)
        
        if (self.selected_iso3 == 'USA') & (us_states):
            self.usstates.plot(edgecolor='k', facecolor='None',
                               alpha=1, linewidth=0.5, zorder=10,
                               ax=ax)
        else:
            # --- Plot Country Border ---
            self.country.plot(edgecolor='k',facecolor='None',
                    alpha=0.7, linewidth=0.5, zorder=11,
                    ax=ax)
        
        ax.axis('off')
        fig.tight_layout()
        
        return fig, ax


class PrettyWorldBaseMap():
    def __init__(self, dpi=250):
        self.countryshapefiles = os.path.join(PATH, 'geometry', 'ne_110m_admin_0_countries', 'ne_110m_admin_0_countries.shp')
        self.worldshapefiles = os.path.join(PATH, 'geometry', 'ne_110m_land', 'ne_110m_land.shp')

        self.dpi = dpi

    def _load_shapes(self):
        self.world = gpd.read_file(self.worldshapefiles)
        self.country = gpd.read_file(self.countryshapefiles)

    def _filter_islands(self, gdf, percentile=0.05):
        gdf["area"] = gdf['geometry'].area / 10**6
        p = gdf["area"].quantile(percentile)
        gdf = gdf.loc[gdf['area'] > p]
        return gdf

    def _filter_centroid(self, gdf, min_lat=-70):
        gdf['centroid'] = gdf['geometry'].centroid
        gdf['centroid'] = [i.coords[0] for i in gdf['centroid']]
        gdf['centroid_lon'] = [i[0] for i in gdf['centroid']]
        gdf['centroid_lat'] = [i[1] for i in gdf['centroid']]

        # --- Filter antarctica ---
        gdf = gdf.loc[gdf['centroid_lat'] > min_lat]

        return gdf

    def _apply_filters(self, gdf):
        gdf = self._filter_islands(gdf)
        gdf = self._filter_centroid(gdf)
        return gdf

    def _filter_shapes(self):
        self.world = self._apply_filters(self.world)
        self.country = self._apply_filters(self.country)

    def fig_ax(self):
        self._load_shapes()
        self._filter_shapes()

        fig, ax = plt.subplots(dpi=self.dpi)
        self.country.plot(edgecolor='k', facecolor='#DADADA',
                          alpha=0.6, linewidth=0.3, ax=ax)
        self.world.plot(edgecolor='#DADADA', facecolor='none',
                        alpha=1, linewidth=0.3, ax=ax)

        ax.axis('off')
        fig.tight_layout()
        return fig, ax
