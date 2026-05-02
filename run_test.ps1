$outputDir = "output"
$exePath = "$outputDir\conjugated_method.exe"

if (!(Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir
}

gcc conjugated_method.c -o $exePath -Wall -Wextra -O2 -lm

if ($LASTEXITCODE -ne 0) {
    Write-Host "Erro na compilacao. O executavel nao foi gerado."
    exit
}

$seed = 42
$tol = "1e-8"
$max_iter = 10000

$bandas = @(7, 27)
$tamanhos = @(1024, 4096, 16384, 65536, 262144, 1048576)

foreach ($banda in $bandas) {
    foreach ($N in $tamanhos) {
        Write-Host "Rodando N=$N banda=$banda"

        & $exePath $N $banda $seed $tol $max_iter

        Write-Host ""
    }
}