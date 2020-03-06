#!/usr/bin/env python
# coding: utf-8

# In[1]:
from threading import Timer
import webbrowser



import sys
sys.path.insert(1, './lib/')
sys.path.insert(1, './test/outputs')
sys.path.insert(1, './test/inputs')
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import base64
import io
import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import json
import numpy as np
import pandas as pd
import wrapper as wp
import datetime 
import astrodynamics  as astro
import datetime as dt
import numpy as np
import jupyterlab_dash
from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory
import os
#
clickprev = 1
nclick_auto_prev = 1
t0        = 0
eclipse   = []
gen_trj   = True
gen_trj_auto = False;
csv = {}
coe = {}
time_plot = []
init_date = 0;
mission_time = 0;
#
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.2)
#
# Default activity profile
#
activity_profile = [{'Name': 'Eclipse', 'T0': 0, 'TF': 100, 'Color': 'yellow'},
                            {'Name': 'Payload', 'T0': 100, 'TF': 200, 'Color': 'red'},
                            {'Name': 'Downlink', 'T0': 200, 'TF': 300, 'Color': 'green'},
                            {'Name': 'Payload', 'T0': 300, 'TF': 400, 'Color': 'red'},
                            {'Name': 'Downlink', 'T0': 400, 'TF': 500, 'Color': 'green'},
                            {'Name': 'Payload', 'T0': 500, 'TF': 600, 'Color': 'red'},
                            {'Name': 'Eclipse', 'T0': 600, 'TF': 700, 'Color': 'yellow'}]
    
#
n = len(activity_profile)
fig.update_xaxes(range=[0, activity_profile[n-1]['TF']], showgrid=True)
fig.update_yaxes(range=[0, 3], row=1, col=1)
#
# Add legends
#
modenames = dict()
#
for i in range(0,n):
      modenames[i] = activity_profile[i]['Name']
#          
for i in range(0,n):
    
         fig.add_trace(go.Scatter(
         x= [(activity_profile[i]['T0'] + activity_profile[i]['TF'])/2.0, (activity_profile[i]['T0'] + activity_profile[i]['TF'])/2.0] ,
         y=[1 , 2],
         hovertemplate = modenames[i],
         showlegend = False,
         mode="text",
         ),row=1, col=1)
    
         fig.add_shape(
         # unfilled Rectangle
         go.layout.Shape(
            type="rect",
            x0=activity_profile[i]['T0'],
            y0=1,
            x1=activity_profile[i]['TF'],
            y1=2,
            line=dict(
                color="Black",
                width=2,
                ),
            fillcolor=activity_profile[i]['Color'],
        ), row=1, col=1)
#       
fig.update_shapes(dict(xref='x', yref='y'))
fig.update_layout(title_text="Activity Profile", title_font_size=30, yaxis_title=" ",xaxis_title= " Elapsed Time (minutes)",title_xanchor='center', title_y= 0.9, title_x=0.5)
#
viewer = jupyterlab_dash.AppViewer()
#
# external JavaScript files
external_scripts = [
    'https://www.google-analytics.com/analytics.js',
    {'src': 'https://cdn.polyfill.io/v2/polyfill.min.js'},
    {
        'src': 'https://cdnjs.cloudflare.com/ajax/libs/lodash.js/4.17.10/lodash.core.js',
        'integrity': 'sha256-Qqd/EfdABZUcAxjOkMi8eGEivtdTkh3b65xCZL4qAQA=',
        'crossorigin': 'anonymous'
    }
]

# external CSS stylesheets
external_stylesheets = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css',
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO',
        'crossorigin': 'anonymous'
    }
]

server = Flask(__name__)
app = dash.Dash(__name__,
                external_scripts=external_scripts,
                external_stylesheets=external_stylesheets,server=server)

app.config['suppress_callback_exceptions'] = True
       
app.layout = html.Div(className='row', children=[
         dcc.Graph(
        id='basic-interactions',
        className='two columns',
        figure=fig,
        config={
            'editable': True},style={'width': '90%', 'display': 'inline-block', 'padding-left':'7%', 'padding-right':'0%'}),
        dcc.Dropdown(
        id='demo-dropdown',
        options=[
            {'label': 'Visibility', 'value': 'visibility'},
            {'label': 'Eclipse', 'value': 'eclipse'},
            {'label': 'RX', 'value': 'rx'},
            {'label': 'RY', 'value': 'ry'},
            {'label': 'RZ', 'value': 'rz'},
            {'label': 'Semi-major Axis', 'value': 'sma'},
            {'label': 'Eccentricity', 'value': 'ecc'},
            {'label': 'Inclination', 'value': 'incl'},
            {'label': 'RAAN', 'value': 'raan'}          
        ],
        value='eclipse',style={'width': '85%', 'padding': '0px 0% 10px 25%'}), 
            html.Div(
            [
        dbc.Button('Generate Activity Profile', id='button',n_clicks=0,style={'width': '200px', 'padding': '0px 0% 0px 0%'}),
        dcc.Upload( id='upload-modes',multiple=False,
        children= ([html.Button('Upload Mode csv', style={'width': '200px', 'padding': '0px 0% 0px 0%'})]))],style={'width': '200px', 'padding': '0px 0% 0px 8%'} ),
                    html.Div(
            [       
        dbc.Button('Generate Trajectory', id='button-trj',n_clicks=0,style={'width': '200px', 'padding': '0px 0% 0px 0%'}),
        dcc.Upload( id='upload-trj',multiple=False,
        children= ([html.Button('Upload Trajectory OEM File', style={'width': '200px', 'padding': '0px 0% 0px 0%'})]))],style={'width': '200px', 'padding': '0px 0% 0px 9%'} ),      
        html.Div(
            [
        dbc.Button('Automatic Generation', id='button-trj-auto',n_clicks=0,style={'width': '200px', 'padding': '0px 0% 0px 0%'}),
        dcc.Input( id='upload-trj-auto', value='Mission_example.json', style={'width': '200px', 'padding': '0px 0% 0px 0%'})],style={'width': '200px', 'padding': '0px 0% 0px 10%'} ),
        html.Div(
            [
        dbc.Button('Generate VTS Files', id='button-vts',n_clicks=0,style={'width': '200px', 'padding': '0px 0% 0px 0%'}), 
        dcc.Input(id='name-vts', value='OUTPUT', type='text',style={'width': '200px', 'padding': '00px 0% 0px 0%'})],style={'width': '200px', 'padding': '0px 0% 0px 11%'} ),      
        html.A(id='download-link',download='SC_Mode.cv',children="Download Mode csv",style={'width': '100000px', 'padding': '0px 10% 91% 75%'}),
        html.Div(id='output-modes-upload'),
        html.Div(
        className='six columns',
        children=[
            html.Div(
                [
                    html.Pre(id='relayout-data'),
                ]
            )
        ]
    )
])

@app.callback(
    Output('basic-interactions', 'figure'),
    [Input('basic-interactions', 'relayoutData'),
     Input('upload-modes', 'contents'), 
     Input('upload-modes', 'filename'),
     Input('button', 'n_clicks'),
     Input('demo-dropdown', 'value'),
     Input('button-trj-auto', 'n_clicks')])
def display_selected_data(relayoutData,contents, names, nclicks, plotvalue,nclicks_auto):
        global fig, n
        global nclick_auto_prev
        global gen_trj_auto
        global modenames
        #
        fig2 = go.Figure(fig) 
        #
        u = 0;
        try : 
            
            fig.update_layout(relayoutData)

            n = len(fig.layout.shapes)
        
            fig.update_shapes(dict(y0=1))
            fig.update_shapes(dict(y1=2))
                 #
                 # Check order of figures
                 #
            xx = np.zeros(n)
            for i in range(0,n):
                    xx[i] =  fig2.layout.shapes[i].x0
                    indd = np.argsort(xx)   

            tend =  fig2.layout.shapes[n-1].x1    
                #
                # Detect if figure have been resized
                #
            pp = 0;
            for i in range(0,n):

                if (not(fig.layout.shapes[indd[i]].x0 == fig2.layout.shapes[indd[i]].x0) and (fig.layout.shapes[indd[i]].x1 == fig2.layout.shapes[indd[i]].x1)) and pp == 0:    
                         if i > 0:
                            fig.layout.shapes[indd[i-1]].x1 = fig.layout.shapes[indd[i]].x0 
                            pp = 1
                            u = 1

                if (not(fig.layout.shapes[indd[i]].x1 == fig2.layout.shapes[indd[i]].x1) and (fig.layout.shapes[indd[i]].x0 == fig2.layout.shapes[indd[i]].x0)) and pp == 0: 
                            if i < n-1:
                                fig.layout.shapes[indd[i+1]].x0 = fig.layout.shapes[indd[i]].x1
                                pp = 1
                                u  = 1
                # 
                # Detect horizontal displacement
                #  
            k = 0; # detect of the shape that has changed has been detected    
            xx = np.zeros(n)
            for i in range(0,n):
                    xx[i] =  fig2.layout.shapes[i].x0

            indd2 = np.argsort(xx) 

            for i in range(0,n):
                if (not(fig.layout.shapes[indd2[i]].x0 == fig2.layout.shapes[indd2[i]].x0) and not(fig.layout.shapes[indd2[i]].x1 == fig2.layout.shapes[indd2[i]].x1)) and pp==0 :    
                        for j in range(0,n-1):
                            if not(indd2[j]==indd2[i]):
                                if ( (fig.layout.shapes[indd2[i]].x0 < fig2.layout.shapes[indd2[j]].x1) and (fig.layout.shapes[indd2[i]].x1 > fig2.layout.shapes[indd2[j+1]].x0)) and pp==0: 
                                     #  
                                    fig.layout.shapes[indd2[i]].x0  = fig2.layout.shapes[indd2[j]].x0
                                    fig.layout.shapes[indd2[i]].x1  = fig2.layout.shapes[indd2[j]].x1
                                    fig.layout.shapes[indd2[i]].y0  = 1
                                    fig.layout.shapes[indd2[i]].y1  = 2
                                    #
                                    fig.layout.shapes[indd2[j]].x0  = fig2.layout.shapes[indd2[i]].x0  
                                    fig.layout.shapes[indd2[j]].x1  = fig2.layout.shapes[indd2[i]].x1 
                                    fig.layout.shapes[indd2[j]].y0  = 1
                                    fig.layout.shapes[indd2[j]].y1  = 2
                                    k = 1 
                                    u = 1
            if k == 0:                            
                    for i in range(0,n):
                        if (not(fig.layout.shapes[indd2[i]].x0 == fig2.layout.shapes[indd2[i]].x0) and not(fig.layout.shapes[indd2[i]].x1 == fig2.layout.shapes[indd2[i]].x1)) and pp==0 :    
                            for j in range(0,n):
                                if not(indd2[j]==indd2[i]):
                                    if ( (fig.layout.shapes[indd2[i]].x0 > fig2.layout.shapes[indd2[j]].x0) and (fig.layout.shapes[indd2[i]].x1 < fig2.layout.shapes[indd2[j]].x1)) and pp==0: 
                                        #
                                        fig.layout.shapes[indd2[i]].x0  = fig2.layout.shapes[indd2[j]].x0
                                        fig.layout.shapes[indd2[i]].x1  = fig2.layout.shapes[indd2[j]].x1
                                        fig.layout.shapes[indd2[i]].y0  = 1
                                        fig.layout.shapes[indd2[i]].y1  = 2
                                        #
                                        fig.layout.shapes[indd2[j]].x0  = fig2.layout.shapes[indd2[i]].x0  
                                        fig.layout.shapes[indd2[j]].x1  = fig2.layout.shapes[indd2[i]].x1 
                                        fig.layout.shapes[indd2[j]].y0  = 1
                                        fig.layout.shapes[indd2[j]].y1  = 2
                                        k = 1 
                                        u = 1
                                        #
            #
            # Check new order and imposed initial and final time limits
            #
            xx = np.zeros(n)
            for i in range(0,n):
                    xx[i] =  fig.layout.shapes[i].x0

            indd = np.argsort(xx) 
            #
            fig.layout.shapes[indd[n-1]].x1 = tend 
            fig.layout.shapes[indd[0]].x0   = 0
            fig.update_shapes(dict(y0=1))
            fig.update_shapes(dict(y0=2))         
            #  
            main_activity_profile(plotvalue)
            #
            for i in range(0,n):
                #
                fig.add_trace(go.Scatter(
                x= [(fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0, (fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0] ,
                y=[1 , 2],
                hovertemplate = modenames[i],
                showlegend = False,
                mode="text"), row=1, col=1)
            #
            fig.update_yaxes(range=[0, 3], row=1, col=1)
            if 'xaxis.range[0]' in relayoutData:
                fig['layout']['xaxis']['range'] = [
                relayoutData['xaxis.range[0]'],
                relayoutData['xaxis.range[1]']
                ]
                fig['layout']['xaxis2']['range'] = [
                relayoutData['xaxis.range[0]'],
                relayoutData['xaxis.range[1]']
                ]
            if 'xaxis2.range[0]' in relayoutData:
                fig['layout']['xaxis']['range'] = [
                relayoutData['xaxis2.range[0]'],
                relayoutData['xaxis2.range[1]']
                ]
                fig['layout']['xaxis2']['range'] = [
                relayoutData['xaxis2.range[0]'],
                relayoutData['xaxis2.range[1]']
                ]
        except:
                fig.update_yaxes(range=[0, 3], row=1, col=1)
                if relayoutData:
                    if 'xaxis.range[0]' in relayoutData:
                        fig['layout']['xaxis']['range'] = [
                        relayoutData['xaxis.range[0]'],
                        relayoutData['xaxis.range[1]']
                        ]
   
                        fig['layout']['xaxis2']['range'] = [
                        relayoutData['xaxis.range[0]'],
                        relayoutData['xaxis.range[1]']
                        ]
                    if 'xaxis2.range[0]' in relayoutData:
                            fig['layout']['xaxis']['range'] = [
                            relayoutData['xaxis2.range[0]'],
                            relayoutData['xaxis2.range[1]']
                            ]
                            fig['layout']['xaxis2']['range'] = [
                            relayoutData['xaxis2.range[0]'],
                            relayoutData['xaxis2.range[1]']
                            ]
            
        if ((contents is not None) and (nclicks is not None)):
            
            activity_profile2 = parse_contents(contents, names, 2)

            global clickprev
            
            if (clickprev==nclicks):
                clickprev = clickprev + 1
                #
                fig = make_subplots(rows=2, cols=1, shared_xaxes=True, shared_yaxes=False,vertical_spacing=0.2) 
                fig.update_shapes(dict(xref='x', yref='y'))
                n   = len(activity_profile2)
                # Set axes properties
                for i in range(0,n):
     
                    fig.add_shape(
                    # unfilled Rectangle
                    go.layout.Shape(
                    type="rect",
                    x0=activity_profile2[i]['T0'],
                    y0=1,
                    x1=activity_profile2[i]['TF'],
                    y1=2,
                    line=dict(
                    color="Black",
                    width=2,
                    ),
                    fillcolor=activity_profile2[i]['Color']
                    ), row=1, col=1)
                    fig.add_trace(go.Scatter(
                    x= [(activity_profile2[i]['T0'] + activity_profile2[i]['TF'])/2.0, (activity_profile2[i]['T0'] + activity_profile2[i]['TF'])/2.0] ,
                    y=[1 , 2],
                    hovertemplate = activity_profile2[i]['Name'],
                    showlegend = False,
                    mode="text"
                    ), row=1, col=1)
                    modenames[i] = activity_profile2[i]['Name']
             
            
            main_activity_profile(plotvalue)
            for i in range(0,n):
            #
                fig.add_trace(go.Scatter(
                x= [(fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0, (fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0] ,
                y=[1 , 2],
                hovertemplate = modenames[i],
                showlegend = False,
                mode="text"), row=1, col=1)
            fig.update_xaxes(range=[0, fig.layout.shapes[n-1].x1], showgrid=True)

            fig.update_layout(title_text="Activity Profile", title_font_size=30, yaxis_title=" ",xaxis_title= " Elapsed Time (minutes)",title_xanchor='center', title_y= 0.9, title_x=0.5)
            fig.update_yaxes(range=[0, 3], row=1, col=1)
        for i in range(0,n):   
                fig.layout.shapes[i].y0 =1
                fig.layout.shapes[i].y1 =2    
        fig.update_yaxes(range=[0, 3], row=1, col=1)
        if relayoutData:
            if 'xaxis.range[0]' in relayoutData:
                fig['layout']['xaxis']['range'] = [
                relayoutData['xaxis.range[0]'],
                relayoutData['xaxis.range[1]']
                ]

                fig['layout']['xaxis2']['range'] = [
                relayoutData['xaxis.range[0]'],
                relayoutData['xaxis.range[1]']
                ]
            if 'xaxis2.range[0]' in relayoutData:
                    fig['layout']['xaxis']['range'] = [
                    relayoutData['xaxis2.range[0]'],
                    relayoutData['xaxis2.range[1]']
                    ]
                    fig['layout']['xaxis2']['range'] = [
                    relayoutData['xaxis2.range[0]'],
                    relayoutData['xaxis2.range[1]']
                    ]
                    
        if (nclicks_auto is not None) and nclicks_auto == nclick_auto_prev:
                #
                nclick_auto_prev = nclick_auto_prev + 1
                main_activity_profile(plotvalue)
                gen_trj = True;
                gen_trj_auto = True;
                
        return fig
    
#
# Function that reads the user-provided spacecraft modes in OEM format
#

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        # Assume that the user uploaded a cvs file
        df = modes_wrapper(io.StringIO(decoded.decode('utf-8')))
        return df
    except Exception as e:
        print(e)
        return html.Div([
            content_string
        ])    
    
def modes_wrapper(input_file):
    global t0
    
    data2 = np.array(pd.read_csv(input_file, header=None, sep=' '))
    
    #f.close()      
    Time   = data2[0:,0] + data2[0:,1]/(3600*24)
    Name   = data2[0:,2:]
    #
    t0     = Time[0]
    Time   = (Time - t0 )*24*60
    t0     = t0*24*60;
    Color  = list()
    nombre = list()
    for j in range(0, len(Time)):
        if  "Payload" in Name[j,0]:
                Color.append('red')
                nombre.append(Name[j,0])
        elif "Downlink" in Name[j,0]:
                Color.append('blue')
                nombre.append(Name[j,0])
        elif "Eclipse" in Name[j,0]:
            Color.append('black')
            nombre.append(Name[j,0])
        else:
            Color.append('yellow')
            nombre.append(Name[j,0])
            
    activity_profile2 = list()
    scmode = dict()
    for j in range(0,len(Time)-1):
        activity_profile2.append({'Name':nombre[j],'T0':Time[j],'TF':Time[j+1],'Color':Color[j]})
    
    return activity_profile2 

def main_activity_profile(Visualization):
    global csv
    global gen_trj
    global time_plot
    global visibility
    global eclipse
    global coe
    global gen_trj_auto
    global modenames
    global n
    global init_date
    #global modenames
    ProjectFile     = "../test/inputs/Mission_example.json"
    TrajectoryFile  = []
    #
    fig.data = ()
    #

    if gen_trj_auto:
        fig.layout.shapes = ()
        time = (time_plot-time_plot[0])/(24*60)
        #
        # Generar secuencia de eventos
        #
        eclipse_change = eclipse[1:-1]- eclipse[0:-2] ;
        #
        eclipse_ind    = np.argwhere(eclipse_change)
        eclipse_ind    = np.concatenate(([0],eclipse_ind[:,0],[len(eclipse)-1]))
        #
        timestart = time[eclipse_ind[0:-1]]
        timeend   = time[eclipse_ind[1:]]
        event     = eclipse[eclipse_ind[1:]]
        eclipse_event = {
        "start" : timestart,
        "end"   : timeend,
        "type"  : event, }
        eclipse_event = pd.DataFrame(eclipse_event) 
        for i in range(0,eclipse_event.shape[0]):
            if eclipse_event.type[i] == 1: 
                dh = eclipse_event.end[i] - eclipse_event.start[i]
                fig.add_shape(
                        # unfilled Rectangle
                        go.layout.Shape(
                        type="rect",
                        x0= (eclipse_event.start[i])*24*60,
                        y0=1,
                        x1=(eclipse_event.start[i]+dh)*24*60,
                        y1=2,
                        line=dict(
                        color="Black",
                        width=2,
                        ),
                        fillcolor='black'
                        ), row=1, col=1)

        #
        # Generar secuencia de visibilidad
        #
        visibility2 = visibility*1*(1-eclipse)
        visibility_change = visibility2[1:-1]- visibility2[0:-2];
        visibility_ind    = np.argwhere(visibility_change)
        visibility_ind    = np.concatenate(([0],visibility_ind[:,0],[len(visibility)-1]))
        timestart = time[visibility_ind[0:-1]]
        timeend   = time[visibility_ind[1:]]
        event     = visibility2[visibility_ind[1:]]*1
        visibility_event = {
        "start" : timestart,
        "end"   : timeend,
        "type"  : event, 
        }
        #
        # Check that there is not eclipse
        #
        visibility_event = pd.DataFrame(visibility_event)
        for i in range(0,visibility_event.shape[0]):
            if visibility_event.type[i] == 1: 
                dh = visibility_event.end[i] - visibility_event.start[i]
                if (visibility_event.start[i] != 1.4):
                    fig.add_shape(
                            # unfilled Rectangle
                            go.layout.Shape(
                            type="rect",
                            x0= (visibility_event.start[i])*24*60,
                            y0=1,
                            x1= (visibility_event.start[i]+dh)*24*60,
                            y1=2,
                            line=dict(
                            color="Black",
                            width=2,
                            ),
                            fillcolor='blue'
                            ), row=1, col=1)

    #
    # Generar secuencia completa de eventos 
    # 
        complete = visibility2 + eclipse;
        complete[complete==2] = 1
        complete_change = complete[1:-1]- complete[0:-2] ;
        complete_ind    = np.argwhere(complete_change)
        complete_ind    = np.concatenate(([0],complete_ind[:,0],[len(complete)-1]))
        timestart = time[complete_ind[0:-1]]
        timeend   = time[complete_ind[1:]]
        event     = complete[complete_ind[1:]]
        complete_event = {
        "start" : timestart,
        "end"   : timeend,
        "type"  : event, 
        }
        complete_event = pd.DataFrame(complete_event)
        for i in range(0,complete_event.shape[0]):
            if complete_event.type[i] == 0: 
                dh = complete_event.end[i] - complete_event.start[i]
                fig.add_shape(
                        # unfilled Rectangle
                        go.layout.Shape(
                        type="rect",
                        x0= (complete_event.start[i])*24*60,
                        y0=1,
                        x1= (complete_event.start[i]+dh)*24*60,
                        y1=2,
                        line=dict(
                        color="Black",
                        width=2,
                        ),
                        fillcolor='red'
                        ), row=1, col=1)
        modenames = dict()
        for i in range(0,n):
            if fig.layout.shapes[i]['fillcolor'] == 'red':
                modenames[i] = 'Payload'
            if fig.layout.shapes[i]['fillcolor'] == 'black':
                modenames[i] = 'Eclipse'
            if fig.layout.shapes[i]['fillcolor']== 'blue':
                modenames[i] = 'Downlink'
            #            
    if gen_trj:
        coe, gs, scmodes, init_date, MissionTime = wp.json_wrapper(ProjectFile)
        t0 = init_date;
        if not TrajectoryFile:
            tf = MissionTime + init_date;
            h  = 10.0;
            csv0 = astro.coe2csv(coe)
            [time,csv] = astro.orbit(csv0,t0,tf,h)
            coe = astro.csv2coe (csv)
            jd0  = init_date
            time = jd0 + time/(3600*24);
        else:
            data = wp.oem_wrapper(TrajectoryFile)
            time = data[0]
            R    = data[1]
            V    = data[2]
            csv = {"RX" : R[:,0],"RY" : R[:,1],"RZ" : R[:,2],"VX" : V[:,0],"VY" : V[:,1], "VZ" : V[:,2]}
            csv = pd.DataFrame(csv) 
            coe = astro.csv2coe(csv)
            #
        eclipse = astro.getElipseCondition(time, csv)
        num_gs = gs.shape[0]
        for i in range(0,num_gs):
            lat    = gs.Latitude[i];
            lon    = gs.Longitude[i]; 
            el     = gs.Elevation[i];        
        visibility = astro.getVisivility(lat, lon, time, csv, el)
        visibility = visibility*1
        time_plot = time*24*60
        time_plot = time_plot-time_plot[0]
    n_plot    = len(Visualization)+1
    plot_curr = Visualization
    if plot_curr == 'rx':
        fig.add_trace(go.Scatter(x=time_plot, y = csv.RX,showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="RX (km)", row=2, col=1)
    elif plot_curr == 'ry':
        fig.add_trace(go.Scatter(x=time_plot, y = csv.RY,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="RY (km)", row=2, col=1)
    elif plot_curr == 'rz':    
        fig.add_trace(go.Scatter(x=time_plot, y = csv.RZ,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="RZ (km)", row=2, col=1)
    elif plot_curr == 'eclipse':  
        fig.add_trace(go.Scatter(x=time_plot, y = eclipse,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="Eclipse", row=2, col=1)
    elif plot_curr == 'visibility':  
        fig.add_trace(go.Scatter(x=time_plot, y = visibility,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="Visibility", row=2, col=1)
    elif plot_curr == 'sma':  
        fig.add_trace(go.Scatter(x=time_plot, y = coe.a,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="Semi-major Axis (Km)", row=2, col=1)
    elif plot_curr == 'ecc':  
        fig.add_trace(go.Scatter(x=time_plot, y = coe.e,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="Eccentricity", row=2, col=1)
    elif plot_curr == 'raan':  
        fig.add_trace(go.Scatter(x=time_plot, y = coe.RA*180/np.pi,
        showlegend = False,marker_color='red'),row=2, col=1),
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="RAAN (deg)", row=2, col=1)
    elif plot_curr == 'incl':  
        fig.add_trace(go.Scatter(x=time_plot, y = coe.incl*180/np.pi,
        showlegend = False,marker_color='red'),row=2, col=1)
        fig.update_xaxes(title_text="Elapsed Time (minutes)", row=2, col=1)
        fig.update_yaxes(title_text="Inclination(deg)", row=2, col=1)
        #
        # 
    n = len(fig.layout.shapes)  
    for i in range(0,n):
            fig.add_trace(go.Scatter(
            x= [(fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0, (fig.layout.shapes[i].x0 + fig.layout.shapes[i].x1)/2.0] ,
            y=[1 , 2],
            hovertemplate = modenames[i],
            showlegend = False,
            mode="text"), row=1, col=1) 
    gen_trj = False            
    gen_trj_auto = False
    return 0

@app.callback(Output('download-link','href'),[Input('button-vts', 'n_clicks'),
    Input('name-vts', 'value')])
def generate_VTS_files(nclicks, filename):
    global n 
    global t0
    global fig
    global init_date
    location = []
    filename2 = []
    t0 = init_date*24*60;
    if (filename is not None) and nclicks >0 : 
            # Order sequence of modes
            xx = np.zeros(n)
            for i in range(0,n):
                    xx[i] =  fig.layout.shapes[i].x0
            indd = np.argsort(xx)
            #
            # Create OEM mode file for VTS
            #
            dateTimeObj = datetime.datetime.now()
            #
            timestampStr = dateTimeObj.strftime("-%d-%b-%Y-%H-%M-%S")
            #
            filename2 = filename + timestampStr +'_OEMVTS.txt'
            path = f"../test/outputs/{filename2}"
                
            with open(path,"w+") as f:
                f.write('CIC_MEM_VERS = 1.0 \n')
                f.write('CREATION_DATE  = 2009-12-08T09:00:00 \n')
                f.write('ORIGINATOR     = SPACEBEL \n')
                f.write('META_START\n')
                f.write('COMMENT Cic infos \n')
                f.write('OBJECT_NAME = CubeSat \n')
                f.write('OBJECT_ID = CubeSat \n')
                f.write('USER_DEFINED_PROTOCOL = CIC \n')
                f.write('USER_DEFINED_CONTENT = SATELLITE_MODES \n')
                f.write('TIME_SYSTEM = UTC \n')
                f.write('META_STOP \n')
                f.write('\n')
                tt = len(fig.layout.shapes)
                for i in range(0,tt):
                        time    = t0 + fig.layout.shapes[indd[i]].x0
                        days    = np.floor(time/(60*24))
                        seconds = (time/(60*24)-days)*3600*24
                        if i == 0:
                            start_day     = days;
                            start_seconds = seconds;
                        
                        f.write(str(int(days))+ ' ' +str(seconds) + ' ' + modenames[indd[i]] + '\n'  ) 
                time    = t0 + fig.layout.shapes[tt-1].x1
                days    = np.floor(time/(60*24))
                seconds = (time/(60*24)-days)*3600*24
                end_day     = days;
                end_seconds = seconds;
                f.write(str(int(days))+ ' ' +str(seconds) + ' ' + "End" + '\n'  )  
            
            #
            f.close() 
            #
            # Create CSV mode file for VTS
            #
            filename3 = filename + timestampStr +'.csv'
            path = f"../test/outputs/{filename3}"
                
            with open(path,"w+") as f:
                tt = len(fig.layout.shapes)
                for i in range(0,tt):
                        time    = t0 + fig.layout.shapes[indd[i]].x0
                        days    = np.floor(time/(60*24))
                        seconds = (time/(60*24)-days)*3600*24
                        if i == 0:
                            start_day     = days;
                            start_seconds = seconds;
                        
                        f.write(str(int(days))+ ' ' +str(seconds) + ' ' + modenames[indd[i]] + '\n'  ) 
                time    = t0 + fig.layout.shapes[tt-1].x1
                days    = np.floor(time/(60*24))
                seconds = (time/(60*24)-days)*3600*24
                end_day     = days;
                end_seconds = seconds;
                f.write(str(int(days))+ ' ' +str(seconds) + ' ' + "End" + '\n'  )  
            

            f.close() 
            
            
            
            filename3 = filename + timestampStr+ 'MAINVTS.vts'
            path2 = f"../test/outputs/{filename3}"
            with open(path2,"w+") as f:
                    f.write('<?xml version="1.0" encoding="UTF-8"?> \n')
                    f.write('<Project Revision="7809">\n')
                    f.write('<General Name="" StartDateTime="')
                    f.write(str(int(start_day))+ ' ' + str(start_seconds) + '" EndDateTime="' + str(int(end_day))+ ' ' + str(end_seconds) + '"/>\n')
                    f.write('<MetaData>\n')
                    f.write(' <Description></Description>\n')
                    f.write('</MetaData>\n')
                    f.write('<MonitorConfiguration>\n')
                    f.write(' <Monitor X="0" Y="0" Height="1030" Width="1920"/>\n')
                    f.write(' <Monitor X="1920" Y="0" Height="1010" Width="1680"/>\n')
                    f.write('</MonitorConfiguration>\n')
                    f.write('<StartOptions TimeRatio="1" UseStateTimeRatio="0" SysTimeSynced="0" Paused="0" Looped="0" Minimized="0" Hidden="0" AutoClosed="0"/>\n')
                    f.write('<TimelineOptions ProjectLocked="1" CursorLocked="0" CursorRatio="0" ViewStart="33282 0.000000" ViewSpan="0" DateFormat="ISODate"/>\n')
                    f.write('<Sky>\n')
                    f.write('<Sun>\n')
                    f.write('  <Prop2d>\n')
                    f.write('   <Icon Anchor="CENTER" Size="MEDIUM" Opacity="100">\n')
                    f.write('    <Font Size="MEDIUM" Color="1 1 1"/>\n')
                    f.write('     <ImageLayer Type="Default"/>\n')
                    f.write('   </Icon>\n')
                    f.write('  </Prop2d>\n')
                    f.write('  <Track Color="0.862745 0.862745 0" PenStyle="SolidLine" PenWidth="2"/>\n')
                    f.write('  <VisibilityCircle ContourColor="0.501961 0.501961 0" FillColor="0 0 0" FillOpacity="50"/>\n')
                    f.write(' </Sun>\n')
                    f.write(' <StarCatalog CatalogMode="Builtin">\n')
                    f.write('  <Track Color="1 1 1" PenStyle="DotLine" PenWidth="1"/>\n')
                    f.write(' </StarCatalog> \n')
                    f.write('</Sky>\n')
                    f.write('<ToBeUsedApps>\n')
                    f.write(' <Application Name="SurfaceView" Id="0" Label="" AutoStarted="1"/>\n')
                    f.write('</ToBeUsedApps>\n')
                    f.write('<Entities>\n')
                    f.write(' <Body Name="Earth" ParentPath="Sol">\n')
                    f.write('  <Prop2d>\n')
                    f.write('   <Icon Anchor="CENTER" Size="MEDIUM" Opacity="100">\n')
                    f.write('    <Font Size="MEDIUM" Color="1 1 1"/>\n')
                    f.write('    <ImageLayer Type="Default"/>\n')
                    f.write('   </Icon>\n')
                    f.write('  </Prop2d>\n')
                    f.write('  <Track Color="0.727337 1 0" PenStyle="SolidLine" PenWidth="2"/>\n')
                    f.write('  <VisibilityCircle ContourColor="0.734829 1 0" FillColor="0.867414 1 0.499992" FillOpacity="60"/>\n')
                    f.write('  <EphemerisMode Mode="Default"/>\n')
                    f.write('  <Layers>\n')
                    f.write('   <BuiltinLayer Name="defaultLayer"/>\n')
                    f.write('  </Layers>\n')
                    f.write(' </Body>\n')
                    f.write('</Entities>\n')
                    f.write('<Events/>\n')
                    f.write('<AdditionalFiles>\n')
                    f.write(' <File Name="')
                    f.write(filename2 + '"/>\n')
                    f.write('</AdditionalFiles>\n')
                    f.write('<States>\n')
                    f.write(' <Instant Time="33282 0" TimeRatio="1" Label="Initial state">\n')
                    f.write('  <AppState Id="0"/>\n')
                    f.write(' </Instant>\n')
                    f.write('</States>\n')
                    f.write('</Project>\n')

            location = urlquote(path)
    return location   


@app.server.route("/test/outputs/<path:path>")
def download(path):
    root_dir = os.getcwd()
    return send_from_directory(
        os.path.join(root_dir, '../test/outputs'), path)

def open_browser(): 
    webbrowser.open_new("http://localhost:{}".format(8898))

if __name__ == '__main__':
           Timer(1, open_browser).start();
           app.run_server(debug=False,port=8898)


# In[2]:





# In[ ]:




