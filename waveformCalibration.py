import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import json
import math
import time
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import os
from decimal import Decimal
from pathlib import Path

configPath = Path(__file__).parent/"config.json"
with open(configPath, 'r') as file:
    config = json.load(file)
defaultFigSize = config['defaultFigSize']
plt.rcParams['lines.linewidth'] = 1

def load_table(filePath, ampCol: list = None):
    LOADERS = {
        ".csv": pd.read_csv,
        ".xlsx": pd.read_excel,
        ".pkl": pd.read_pickle,
        ".p": pd.read_pickle
    }

    ext = Path(filePath).suffix.lower()

    if ext not in LOADERS:
        raise ValueError(
            f"Unsupported file format '{ext}'. "
            f"Supported formats: {', '.join(LOADERS.keys())}"
        )

    # Load data
    data = LOADERS[ext](filePath)

    # Sort by wavelength column (index 0)
    data = data.sort_values(by=data.columns[0])

    # Determine which amplitude columns to load
    if ampCol is None:
        amp_indices = list(range(1, data.shape[1]))  # all except wl column
    else:
        amp_indices = ampCol

    # Build output DataFrame (wl + selected amp columns)
    selected_cols = [0] + amp_indices
    data_selected = data.iloc[:, selected_cols]

    # Extract channel names (exclude wl column)
    try:
        channelNames = [str(data.columns[i]) for i in amp_indices]
    except Exception:
        channelNames = None

    # Build name structure
    name = [Path(filePath).stem, channelNames]

    return data_selected, name

def selectFile(*fileTypes, initPath=None, titleText="Select a file"):
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)

    # 1. Handle Path Logic
    i_dir = None
    i_file = None

    if initPath:
        p = Path(initPath).resolve()
        if p.is_dir():
            i_dir = str(p)
        else:
            i_dir = str(p.parent)
            i_file = p.name
            
    # 2. Build Extension Library
    extensionLib = {
        '.p': ("Pickle files", "*.p *.pickle"),
        '.xlsx': ("Excel Workbook", "*.xlsx"),
        '.csv': ("CSV files", "*.csv"),
        '.json': ("JSON files", "*.json"),
    }

    selected_types = []
    if fileTypes:
        for ext in fileTypes:
            ext = ext.lower()
            if ext in extensionLib:
                selected_types.append(extensionLib[ext])
            else:
                selected_types.append((f"{ext.upper()} files", f"*{ext}"))
    else:
        selected_types = []
    
    selected_types.append(("All files", "*.*"))

    # 3. Open Dialog
    filePath = filedialog.askopenfilename(
        title=titleText,
        filetypes=selected_types,
        initialdir=i_dir,
        initialfile=i_file
    )

    root.destroy()
    return filePath

def calFitScore(dataIntensity, refIntensity, option='SSD'):
    axis=0
    if max( len(dataIntensity.shape), len(refIntensity.shape) ) ==2:
        axis=1
    if option=='sum difference':
        score = -np.sum(np.abs(refIntensity-dataIntensity), axis=axis)
    elif option=='SSD': # sum of squared differences
        score = -np.sum((refIntensity-dataIntensity)**2, axis=axis)
    elif option=='CCF': # cross-correlation function
        score = np.sum(dataIntensity*refIntensity, axis=axis)
    else:
        raise ValueError
    return score

def select_data(domain, intensity, selectionRange, maxShift=0, option="extend"): # include one point left of the range, and another point to the right
    # option: "extend", "strict"
    start = selectionRange[0]-maxShift
    stop  = selectionRange[1]+maxShift

    startIndex = -1
    stopIndex = len(domain)+1
    for i in range(len(domain)):
        if domain[i]<start:
            startIndex = i
        if domain[i]>stop:
            stopIndex = i
            break
    if option=="strict":
        startIndex += 1
        stopIndex  -= 1
    if startIndex==-1:
        startIndex=0
    if stopIndex==len(domain)+1:
        stopIndex=len(domain)
    return domain[startIndex: stopIndex], intensity[startIndex: stopIndex] # returning pointer instead of copied array, should actually be copied array

def findShift(scanWavelength, scanData, refWavelength, refIntensity, bandRange):
    '''
    This function slides observed spectrum over the reference waveform.
    - risk of upsampling
    - faster than findShift2
    '''
    shiftResolution = config['shiftResolution']
    maxShift        = config['maxShift']

    decimalPoint = len(str(Decimal(str(shiftResolution))).split('.')[1])
    startRange   = bandRange[0] - maxShift*2
    stopRange    = bandRange[1] + maxShift*2
    startRange   = round(math.ceil(startRange/shiftResolution)*shiftResolution, decimalPoint)
    stopRange    = round(math.floor(stopRange/shiftResolution)*shiftResolution, decimalPoint)

    n = round( (stopRange-startRange)/shiftResolution+1 )
    grid = np.linspace(startRange, stopRange, n, endpoint=True) # grid is as large as bandRange +- maxShift

    fref = interp1d(refWavelength, refIntensity, kind='quadratic')
    refIntensity = fref(grid)
    refWavelength, refIntensity = select_data(grid, refIntensity, bandRange, 0, option='strict')

    fscan = interp1d(scanWavelength, scanData, kind='quadratic')
    scanData = fscan(grid)

    selectWavelength, selectData = select_data(grid, scanData, bandRange, 0, option='strict')
    startIndex = np.where(grid == selectWavelength[0])[0][0]
    stopIndex  = np.where(grid == selectWavelength[-1])[0][0]
    nShift     = round( maxShift//shiftResolution )

    shiftIndex = np.linspace(-nShift, nShift, 2*nShift+1, endpoint=True, dtype=np.int16)
    shiftIntensities = np.full((len(shiftIndex), stopIndex-startIndex+1), np.nan)

    for i in range(len(shiftIndex)):
        index        = -shiftIndex[i] # The intensities shift in the opposite side relative to the wavelength domain.
        selectedData = scanData[startIndex+index: stopIndex+index+1]
        shiftIntensities[i] = selectedData
    scoreArray   = calFitScore(shiftIntensities, refIntensity, option='SSD')
    maxScore     = np.max(scoreArray)
    optimumIndex = np.where(scoreArray==maxScore)[0][0]
    optimumShift = shiftIndex[optimumIndex]*shiftResolution
    
    scoreDataSet = pd.DataFrame()
    scoreDataSet['shiftValue'] = shiftIndex*shiftResolution
    scoreDataSet['score']      = scoreArray
    return scoreDataSet, optimumShift

def findShift2(scanWavelength, scanData, refWavelength, refIntensity, bandRange):
    '''
    This function slides reference waveform over observed spectrum.
    '''
    shiftResolution = config['shiftResolution']
    maxShift        = config['maxShift']
    refWavelength , refIntensity = select_data( refWavelength, refIntensity, bandRange, maxShift=maxShift, option="extend")
    scanWavelength,     scanData = select_data(scanWavelength,     scanData, bandRange, maxShift=maxShift, option="strict")

    fref = interp1d(refWavelength, refIntensity, kind='cubic')

    nShift      = round( maxShift//shiftResolution )
    shiftIndices = np.linspace(-nShift, nShift, 2*nShift+1, endpoint=True, dtype=np.int16)
    shiftValues = shiftIndices*shiftResolution

    scoreArray = []
    for val in shiftValues:
        shifted_wl = scanWavelength + val
        obs_points, obs_intensity = select_data(shifted_wl, scanData, bandRange, maxShift=0, option='strict')
        nPoints = len(obs_points)
        ref_interpolatedIntensity = fref(obs_points)

        score = calFitScore(obs_intensity, ref_interpolatedIntensity, option='CCF')
        scoreArray.append(score/nPoints)
    maxScore     = np.max(scoreArray)
    optimumIndex = np.where(scoreArray==maxScore)[0][0]
    optimumShift = shiftValues[optimumIndex]
    
    scoreDataSet = pd.DataFrame()
    scoreDataSet['shiftValue'] = shiftValues
    scoreDataSet['score']      = scoreArray
    return scoreDataSet, optimumShift

def plotOriginal(scanWavelength, scanData, refWavelength, refIntensity):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=defaultFigSize)
    
    # First plot
    bandRange1 = config['bandRange1']
    domain1, data1 = select_data(scanWavelength,     scanData, bandRange1, 0)
    refDom1, ref1  = select_data( refWavelength, refIntensity, bandRange1, 0)
    ax1.plot(domain1, data1, label = "observed spectrum")
    ax1.plot(refDom1,  ref1, label = "reference band 1")
    ax1.legend(loc="lower right")

    # Second plot
    bandRange2 = config['bandRange2']
    domain2, data2 = select_data(scanWavelength, scanData, bandRange2,0)
    refDom2,  ref2 = select_data( refWavelength, refIntensity, bandRange2,0)
    ax2.plot(domain2, data2, label = "observed spectrum")
    ax2.plot(refDom2, ref2, label = "reference band 2", color="limegreen")
    ax2.legend(loc="lower right")

    plt.tight_layout()
    plt.show()
    return None

def plotScores(scoreData1, scoreData2):
    shiftValue1 = scoreData1['shiftValue']
    score1      = scoreData1['score']
    shiftValue2 = scoreData2['shiftValue']
    score2      = scoreData2['score']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=defaultFigSize)
    circleSize = 50

    ax1.plot(shiftValue1, score1, label='band 1')
    maxScore1 = np.max(score1)
    maxIndex_1 = np.where(score1==maxScore1)[0][0]
    ax1.scatter(shiftValue1[maxIndex_1], maxScore1, s=circleSize, facecolors='none', edgecolors='red', linewidths=2, marker='o')
    ax1.set_xlabel('lags (nm)')
    ax1.legend(loc="lower right")

    ax2.plot(shiftValue2, score2, label='band 2')
    maxScore2  = np.max(score2)
    maxIndex_2 = np.where(score2==maxScore2)[0][0]
    ax2.scatter(shiftValue2[maxIndex_2], maxScore2, s=circleSize, facecolors='none', edgecolors='red', linewidths=2, marker='o')
    ax2.set_xlabel('lags (nm)')
    ax2.legend(loc="lower right")

    plt.tight_layout()
    plt.show()
    return None

def plotCorrectedWaveform(scanWavelength, scanData, refWavelength, refIntensity):
    plt.figure(figsize=defaultFigSize)
    plt.plot(scanWavelength, scanData, label='corrected waveform')
    plt.plot(refWavelength, refIntensity, label='reference')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.title('Corrected waveform with reference')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.grid(True)
    plt.show()

def plotBandShifts(obs_wl_shifted, obs_flux, ref_wl, ref_flux, shift1, shift2):
    range1 = config["bandRange1"]
    range2 = config["bandRange2"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    ax1.plot(ref_wl, ref_flux, label='reference')
    ax1.plot(obs_wl_shifted.copy()+shift1, obs_flux, label='B-band')
    ax1.set_ylabel("I (A.U.)")
    ax1.legend(loc=1)
    ax1.set_xlim(*range1)
    ax1.grid(True)

    ax2.plot(ref_wl, ref_flux, label='reference')
    ax2.plot(obs_wl_shifted.copy()+shift2, obs_flux, label='B-band')
    ax2.set_xlabel("wavelength (nm)")
    ax2.set_ylabel("I (A.U.)")
    ax2.legend(loc=1)
    ax2.set_xlim(*range2)
    ax2.grid(True)

    fig.suptitle(f"Band shifts applied separately")
    plt.tight_layout()
    plt.show()
    pass

def findCoefficients(scanWavelength, scanData, refWavelength, refIntensity, print_to_terminal=True):
    bandRange1 = config['bandRange1']
    bandRange2 = config['bandRange2']
    if config['showPlots']:
        plotOriginal(scanWavelength, scanData, refWavelength, refIntensity)
        pass

    score1, shift1 = findShift2(scanWavelength, scanData, refWavelength, refIntensity, bandRange1)
    score2, shift2 = findShift2(scanWavelength, scanData, refWavelength, refIntensity, bandRange2)
    if print_to_terminal:
        print(f'shift 1 = {float(shift1)}')
        print(f'shift 2 = {float(shift2)}')
    if config['showPlots']:
        plotBandShifts(scanWavelength, scanData, refWavelength, refIntensity, shift1, shift2)
        pass

    midRange1 = sum(bandRange1)/2
    midRange2 = sum(bandRange2)/2
    k = float((shift2-shift1)/(midRange2-midRange1))
    b = float(shift1 - k*midRange1)

    if print_to_terminal:
        print(f'{k= }')
        print(f'{b= }')

    newWavelength = scanWavelength + scanWavelength*k+b
    if config['showPlots']:
        plotScores(score1, score2)
        plotCorrectedWaveform(newWavelength, scanData, refWavelength, refIntensity)
    return newWavelength, (k, b)

def cal_avg(wavelengths, spectra):
    min_wl = wavelengths.min().max()
    max_wl = wavelengths.max().min()
    selectionRange = [min_wl, max_wl]

    nPoints = []
    for dataset in wavelengths.columns:
        wl, intensity = select_data(wavelengths[dataset], spectra[dataset], selectionRange, maxShift=0, option="strict")
        nPoints.append( (len(wl), str(dataset)) )
    nPoints.sort()
    templateName = nPoints[0][1]

    templateGrid, intensity = select_data(wavelengths[templateName], spectra[templateName], selectionRange, maxShift=0, option="strict")
    
    spectraHolder = pd.DataFrame()
    for scan in spectra:
        fscan = interp1d(wavelengths[scan], spectra[scan], kind='cubic')
        spectraHolder[scan] = fscan(templateGrid)
    array = spectraHolder.to_numpy()
    avg_spectrum = np.mean(array, axis=1)
    avg_corrected_waveform = pd.DataFrame()
    avg_corrected_waveform['Wavelength_(nm)'] = templateGrid
    avg_corrected_waveform['avg_intensity']   = avg_spectrum
    return avg_corrected_waveform

def format_elapsed(seconds):
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        m, s = divmod(seconds, 60)
        return f"{int(m)} min {s:.1f} sec"
    else:
        h, rem = divmod(seconds, 3600)
        m, s = divmod(rem, 60)
        return f"{int(h)} hr {int(m)} min {s:.1f} sec"

def selectDataset(columns):
    columns = list(columns)
    print(f'available data set: {', '.join(columns)}')
    print('Dataset selection')
    print('  - if select all, leave blank.')
    print("""  - use comma ',' to separate column name""")

    includes = input(" : ")
    if not len(includes.strip()):
        print('All columns are selected.')
        return columns
    includes = includes.split(',')
    includes = [e.strip().strip('''"''').strip("""'""").strip() for e in includes]

    selectedColumns = [x for x in columns if x in includes]
    print(f'selected: {', '.join(selectedColumns)}')

    return selectedColumns

def main(showPlot=False, outputDir=None):
    dataFileName = selectFile(initPath=Path(__file__).parent/'uncalibrated_spectrum', titleText="Select observed spectrum data file")
    refFileName  = selectFile(initPath=Path(__file__).parent/'reference_spectrum', titleText="Select reference data file")
    startTime     = time.time()
    print("")
    print(f'observed data file: {dataFileName}')
    print(f'reference file    : {refFileName}')
    spectrumTable, spectrumName = load_table(dataFileName)
    refFile, _ = load_table(refFileName)

    spectrumTable = spectrumTable.dropna(axis=1, how='all')
    spectrumTable = spectrumTable.drop_duplicates()
    spectrumTable = spectrumTable.sort_values(by=spectrumTable.columns[0])
    spectrumTable = spectrumTable.reset_index(drop=True)

    refFile       = refFile.dropna(axis=1, how='all')
    refFile       = refFile.drop_duplicates()
    refFile       = refFile.sort_values(by=refFile.columns[0])
    refFile       = refFile.reset_index(drop=True)
    refWavelength = refFile.iloc[:, 0].to_numpy()
    refIntensity  = refFile.iloc[:, 1].to_numpy()

    I_scale = np.mean(spectrumTable.iloc[:,1]) / np.mean(refIntensity)
    refIntensity *= I_scale

    selectedSignal  = selectDataset(spectrumTable.columns[1:])

    k_list = []
    b_list = []
    newWavelengths = pd.DataFrame()
    for scanName in selectedSignal:
        print("\ncalculating: "+scanName)
        obsWavelength = spectrumTable.iloc[:, 0].to_numpy(copy=True)
        currentSpectra = spectrumTable[scanName].to_numpy(copy=True)
        calibratedWavelength, coefficients = findCoefficients(obsWavelength, currentSpectra, refWavelength, refIntensity)
        newWavelengths[scanName] = calibratedWavelength
        k, b = coefficients
        k_list.append(k)
        b_list.append(b)
    
    columnList = spectrumTable.columns
    for column in columnList:
        if column not in selectedSignal:
            spectrumTable = spectrumTable.drop(column, axis=1)

    avg_corrected_waveform = cal_avg(newWavelengths, spectrumTable)
    
    if config['saveCorrectedSignal']:
        result_dir = Path(outputDir)/spectrumName[0]
        result_dir.mkdir(parents=True, exist_ok=True)
        print("\nsaving files")
        spectrumFileName = result_dir / "corrected_intensities.csv"
        wavelengthsFileName = result_dir / "corrected_wavelengths.csv"
        coefFileName = result_dir / 'correction_coefficients.csv'
        avgFileName = result_dir / "avg_corrected_waveform.csv" if len(selectedSignal)>1 else result_dir / f'{selectedSignal[0]}_corrected_waveform.csv'
        spectrumTable.to_csv(spectrumFileName, index=False, header=True)
        print(f' - {spectrumFileName}')
        newWavelengths.to_csv(wavelengthsFileName, index=False, header=True)
        print(f' - {wavelengthsFileName}')
        correction_coefficients = pd.DataFrame()
        correction_coefficients['signalName']=selectedSignal
        correction_coefficients['slope'] = k_list
        correction_coefficients['offset'] = b_list
        correction_coefficients.to_csv(coefFileName, index=False, header=True)
        print(f' - {coefFileName}')
        avg_corrected_waveform.to_csv(avgFileName, index=False, header=True)
        print(f' - {avgFileName}')

    finishTime = time.time()
    elapsed = finishTime - startTime
    print("\nElapsed time:", format_elapsed(elapsed))

    if showPlot:  # final result
        plt.plot(avg_corrected_waveform.iloc[:, 0], avg_corrected_waveform.iloc[:, 1], label='averaged corrected spectrum')
        plt.plot(refWavelength, refIntensity, label='reference spectrum')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity')
        plt.title(f'Calibrated {spectrumName[0]} Waveform')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return None

if __name__=="__main__":
    print('\nInitializing...')
    outDir = Path(__file__).parent/'calibrated_output_spectrum'
    main(showPlot=False, outputDir=outDir)