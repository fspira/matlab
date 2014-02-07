xCR = [10 18 16 14]
yCR = [1.36 1.20 1.15 1.49]
plot (xCR(:), yCR(:),'x')
plot (xCR(:), yCR(:),'x')
axis([0 25 0.8 1.8])
 xlabel ('Contractile ring diameter [µm]','FontSize', 16);
            ylabel('Ratio P-S [A.U.]','FontSize', 16);
            title('Fluorescence ratio relative to cleavage furrow diameter','FontSize', 16);