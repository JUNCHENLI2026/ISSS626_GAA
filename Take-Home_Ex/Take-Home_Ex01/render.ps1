param([switch]$Verify)
$ErrorActionPreference = 'Stop'

# R on Windows needs a Windows-supported UTF-8 locale for non-ASCII paths.
# These settings apply only while this script runs.
$previousLocation = Get-Location
$environmentNames = @('LC_ALL', 'LC_CTYPE', 'LANG', 'TEMP', 'TMP')
$previousEnvironment = @{}
foreach ($name in $environmentNames) {
    $previousEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
}

try {
    Set-Location -LiteralPath $PSScriptRoot
    $env:LC_ALL = 'English_United States.utf8'
    $env:LC_CTYPE = 'English_United States.utf8'
    $env:LANG = 'en_US.UTF-8'
    $temporaryDirectory = Join-Path $PSScriptRoot '.render-tmp'
    New-Item -ItemType Directory -Force -Path $temporaryDirectory | Out-Null
    $env:TEMP = $temporaryDirectory
    $env:TMP = $temporaryDirectory

    if ($Verify) {
        $rCommand = Get-Command Rscript -ErrorAction SilentlyContinue
        if ($rCommand) {
            $rExecutable = $rCommand.Source
        } else {
            $rInstallation = Get-ItemProperty 'HKLM:\SOFTWARE\R-core\R' -ErrorAction Stop
            $rExecutable = Join-Path $rInstallation.InstallPath 'bin\Rscript.exe'
        }
        foreach ($check in @('verify_revision.R', 'verify_robustness.R')) {
            & $rExecutable (Join-Path 'data/Ketapang' $check) 'data/Ketapang'
            if ($LASTEXITCODE -ne 0) { throw "Verification failed: $check" }
        }
    }

    $quartoExecutable = (Get-Command quarto -ErrorAction Stop).Source
    foreach ($document in @('Take-Home_Ex01.qmd', 'executive_summary.qmd', 'index.qmd')) {
        & $quartoExecutable render $document
        if ($LASTEXITCODE -ne 0) { throw "Rendering failed: $document" }
    }
} finally {
    Set-Location -LiteralPath $previousLocation.Path
    foreach ($name in $environmentNames) {
        [Environment]::SetEnvironmentVariable($name, $previousEnvironment[$name], 'Process')
    }
}
