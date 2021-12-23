
$files = 
    "RR.ps1",
    "yakpackage\manifest.yml",
    "MinSurfaceGH/MinSurfacev7GH.csproj",
    # "MinSurfaceCommands\MinSurfaceCommands.vcxproj",
    "MinSurfaceGH/Properties/AssemblyInfo.cs";


ForEach($file in $files) {
    if (!(Test-Path $file)) {
        throw "File " + $file + " not found!";
    }
    Write-Host $file -ForegroundColor Green;
    (Get-Content $file) | Select-String -Pattern "0.1" | Write;
}
       
