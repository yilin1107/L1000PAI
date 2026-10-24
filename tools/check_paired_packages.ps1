[CmdletBinding()]
param(
    [string]$CorePath = '',
    [string]$DataPath = '',
    [string]$RExecutable = '',
    [switch]$Quick
)

$ErrorActionPreference = 'Stop'

function Resolve-RExecutable {
    param([string]$Requested)

    if ($Requested) {
        return (Resolve-Path -LiteralPath $Requested).Path
    }

    $fromPath = Get-Command R.exe -ErrorAction SilentlyContinue
    if ($fromPath) {
        return $fromPath.Source
    }

    $candidates = Get-ChildItem -Path 'C:\Program Files\R' -Recurse -Filter R.exe -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -match '[\\/]bin[\\/]R\.exe$' } |
        Sort-Object FullName -Descending
    if ($candidates) {
        return $candidates[0].FullName
    }

    throw 'R.exe was not found. Pass -RExecutable with an explicit path.'
}

function Copy-PackageAllowlist {
    param(
        [string]$Source,
        [string]$Destination
    )

    $entries = @(
        '.Rbuildignore', 'DESCRIPTION', 'NAMESPACE', 'LICENSE', 'LICENSE.md', 'NEWS.md',
        'README.md', 'R', 'man', 'inst', 'tests', 'vignettes'
    )
    New-Item -ItemType Directory -Force -Path $Destination | Out-Null
    foreach ($entry in $entries) {
        $sourceEntry = Join-Path $Source $entry
        if (Test-Path -LiteralPath $sourceEntry) {
            Copy-Item -LiteralPath $sourceEntry -Destination $Destination -Recurse -Force
        }
    }
}

function Invoke-RCommand {
    param(
        [string]$Executable,
        [string[]]$Arguments
    )

    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "R command failed with exit code ${LASTEXITCODE}: $($Arguments -join ' ')"
    }
}

if (-not $CorePath) {
    $CorePath = Join-Path $PSScriptRoot '..'
}
$core = (Resolve-Path -LiteralPath $CorePath).Path
if (-not $DataPath) {
    if ($env:L1000PAI_DATA_DIR) {
        $DataPath = $env:L1000PAI_DATA_DIR
    } else {
        $DataPath = Join-Path (Split-Path $core -Parent) 'L1000PAIdata'
    }
}
$data = (Resolve-Path -LiteralPath $DataPath).Path
$rExe = Resolve-RExecutable -Requested $RExecutable
$rScript = Join-Path (Split-Path $rExe -Parent) 'Rscript.exe'
if (-not (Test-Path -LiteralPath $rScript)) {
    $rScript = [IO.Path]::ChangeExtension($rExe, $null) + 'script.exe'
}
if (-not (Test-Path -LiteralPath $rScript)) {
    throw "Rscript.exe was not found beside $rExe"
}

$auditRoot = Join-Path ([IO.Path]::GetTempPath()) ("L1000PAI-paired-check-" + [guid]::NewGuid().ToString('N'))
$coreStage = Join-Path $auditRoot 'src\L1000PAI'
$dataStage = Join-Path $auditRoot 'src\L1000PAIdata'
$isolatedLibrary = Join-Path $auditRoot 'library'
$dataCheckWork = Join-Path $auditRoot 'data-check'
$checkWork = Join-Path $auditRoot 'check'
New-Item -ItemType Directory -Force -Path $isolatedLibrary, $dataCheckWork, $checkWork | Out-Null

Copy-PackageAllowlist -Source $core -Destination $coreStage
Copy-PackageAllowlist -Source $data -Destination $dataStage

$existingLibraries = & $rScript -e 'cat(.libPaths(), sep=.Platform$path.sep)'
if ($LASTEXITCODE -ne 0) {
    throw 'Unable to query the active R library paths.'
}
$pathSeparator = [IO.Path]::PathSeparator
$previousRLibs = $env:R_LIBS
$previousAuditLibrary = $env:L1000PAI_AUDIT_LIBRARY
$env:R_LIBS = $isolatedLibrary + $pathSeparator + $existingLibraries
$env:L1000PAI_AUDIT_LIBRARY = $isolatedLibrary

try {
    Push-Location $auditRoot

    $dataBuild = @('CMD', 'build')
    if ($Quick) { $dataBuild += '--no-build-vignettes' }
    $dataBuild += $dataStage
    Invoke-RCommand -Executable $rExe -Arguments $dataBuild

    $dataTarball = Get-ChildItem -LiteralPath $auditRoot -Filter 'L1000PAIdata_*.tar.gz' |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $dataTarball) { throw 'The L1000PAIdata source tarball was not created.' }

    Set-Location $dataCheckWork
    $dataCheckArguments = @('CMD', 'check')
    if ($Quick) { $dataCheckArguments += @('--no-manual', '--ignore-vignettes') }
    $dataCheckArguments += $dataTarball.FullName
    Invoke-RCommand -Executable $rExe -Arguments $dataCheckArguments

    Set-Location $auditRoot
    Invoke-RCommand -Executable $rExe -Arguments @(
        'CMD', 'INSTALL', "--library=$isolatedLibrary", $dataTarball.FullName
    )

    Invoke-RCommand -Executable $rScript -Arguments @(
        '-e',
        "installed <- normalizePath(find.package('L1000PAIdata'), winslash='/'); expected <- normalizePath(Sys.getenv('L1000PAI_AUDIT_LIBRARY'), winslash='/'); stopifnot(startsWith(installed, expected), packageVersion('L1000PAIdata') >= '1.0.0'); v <- L1000PAIdata::verify_l1000pai_data_files(stop_on_error = TRUE); cat('verified cached libraries:', sum(v[['exists']]), '\n')"
    )

    $coreBuild = @('CMD', 'build')
    if ($Quick) { $coreBuild += '--no-build-vignettes' }
    $coreBuild += $coreStage
    Invoke-RCommand -Executable $rExe -Arguments $coreBuild

    $coreTarball = Get-ChildItem -LiteralPath $auditRoot -Filter 'L1000PAI_*.tar.gz' |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $coreTarball) { throw 'The L1000PAI source tarball was not created.' }

    Set-Location $checkWork
    $checkArguments = @('CMD', 'check')
    if ($Quick) { $checkArguments += @('--no-manual', '--ignore-vignettes') }
    $checkArguments += $coreTarball.FullName
    Invoke-RCommand -Executable $rExe -Arguments $checkArguments

    & $rScript -e "writeLines(capture.output(sessionInfo()), '$($auditRoot.Replace('\', '/'))/sessionInfo.txt')"
    if ($LASTEXITCODE -ne 0) { throw 'Unable to write sessionInfo.txt.' }
    Get-FileHash -Algorithm SHA256 $dataTarball.FullName, $coreTarball.FullName |
        Format-Table -AutoSize |
        Out-String |
        Set-Content -LiteralPath (Join-Path $auditRoot 'artifact_hashes.txt') -Encoding UTF8

    Write-Host "Paired package check completed. Audit directory: $auditRoot"
} finally {
    Pop-Location -ErrorAction SilentlyContinue
    $env:R_LIBS = $previousRLibs
    $env:L1000PAI_AUDIT_LIBRARY = $previousAuditLibrary
}
