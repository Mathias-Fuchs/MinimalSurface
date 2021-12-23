# my "CI workflow":
# 1) use printOccuringVersions.ps1 to check all version numbers are correct
# 2) use RR.ps1 to build the release
# 3) execute this script in dev powershell 

Yak.exe search MinSurface | sed 's/.*(//g' | sed 's/).*//g' > publishedversion.txt
$publishee = "yakpackage/minsurface-" + $args[0] + "-rh7_3-win.yak"
if (!(Test-Path $publishee)) {throw "Publishee doesn't exist";}
Yak.exe push $publishee

"$args[0]" | Out-File publishedversion.txt

.\rename.ps1 $args[0] $args[1]

New-Item -Path ..\..\..\AppData\Roaming\McNeel\Rhinoceros\packages\7.0\MinSurface -Name $args[1] -ItemType "directory";

.\printOccuringVersions.ps1
