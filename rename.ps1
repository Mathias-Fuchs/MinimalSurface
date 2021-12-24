
$files = 
 "RR.ps1",
 "yakpackage\manifest.yml",
 "MinSurfaceGH\MinSurfacev7GH.csproj",
 "MinSurfaceCommands\MinSurfaceCommands.vcxproj",
 "MinSurfaceCommands\MinSurfaceCommands.rc",
 "MinSurfaceGH\Properties\AssemblyInfo.cs";

ForEach($file in $files) {
    if (!(Test-Path $file)) {throw "File " + $file + " not found!";}
    Write-Host $file -ForegroundColor Green;
    (Get-Content $file).Replace($args[0], $args[1]) | Set-Content $file
}
       
