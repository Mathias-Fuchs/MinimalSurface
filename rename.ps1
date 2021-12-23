
$files = "RR.ps1", "yakpackage\manifest.yml", "MinSurfaceGH/MinSurfacev7GH.csproj";

            
ForEach($file in $files) {
    Test-Path $file
    (Get-Content $file).replace($args[0], $args[1]) | Set-Content $file
}
       
