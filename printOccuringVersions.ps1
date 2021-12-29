
$files = 
    "RR.ps1",
    "yakpackage\manifest.yml",
    "MinSurfaceGH\MinSurfacev7GH.csproj",
    "MinSurfaceCommands\MinSurfaceCommands.vcxproj",
    "MinSurfaceCommands\MinSurfaceCommandsPlugIn.cpp",
    "MinSurfaceCommands\MinSurfaceCommands.rc",
    "MinSurfaceGH\Properties\AssemblyInfo.cs";

ForEach($file in $files) {
    if (!(Test-Path $file)) {
        throw "File " + $file + " not found!";
    }
    Write-Host $file -ForegroundColor Green;
    (Get-Content $file) | Select-String -Pattern "0.1" | Write;
}

Write-Host "Gha assembly:" -ForegroundColor Green;
(Get-Item yakpackage\*.gha).VersionInfo | Write-Output;

Write-Host "Native dll:" -ForegroundColor Green;
dumpbin /dependents yakpackage\*.rhp
dumpbin /headers yakpackage\*.rhp
