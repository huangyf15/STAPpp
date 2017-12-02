fid1=fopen("rect002.txt","wt");
strings_o=["xdir[0]" "xdir[1]" "xdir[2]" "ydir[0]" "ydir[1]" "ydir[2]" "zdir[0]" "zdir[1]" "zdir[2]"];
for j=1:24
    for i=j:-1:1
        modi=mod(i,6);
        modj=mod(j,6);
        if (modi==0)
            modi = 6;
        end
        if (modj ==0)
            modj = 6;
        end
        parti=fix(i/6); %begin from 0
        partj=fix(j/6);
        if (modi>=1 && modi<=3)
            if (modj>=1 && modj<=3)
                fprintf(fid1,'%s%i%s%i%s%s%s%s%s\n',"Matrix[",j*(j+1)/2-i,"] = MatrixLit[",vecpos(3*parti+1,3*partj+1),"] * ",strings_o(6+modi),"*",strings_o(6+modj),";");
            else
                fprintf(fid1,'%s%i%s%i%s%s%s%s%s%i%s%s%s%s%s\n',"Matrix[",j*(j+1)/2-i,"] = MatrixLit[",vecpos(3*parti+1,3*partj+2),"] * ",strings_o(6+modi),"*",strings_o(modj-3)," + MatrixLit[",vecpos(3*parti+1,3*parti+3),"] * ",strings_o(6+modi),"*",strings_o(modj),";");
            end
        else
            if (modj>=1 && modj <=3)
                fprintf(fid1,'%s%i%s%i%s%s%s%s%s%i%s%s%s%s%s\n',"Matrix[",j*(j+1)/2-i,"] = MatrixLit[",vecpos(3*parti+2,3*partj+1),"] * ",strings_o(modi-3),"*",strings_o(6+modj)," + MatrixLit[",vecpos(3*parti+3,3*partj+1),"] * ",strings_o(modi),"*",strings_o(6+modj),";");
            else
                fprintf(fid1,'%s%i%s%i%s%s%s%s%s',"Matrix[",j*(j+1)/2-i,"] = MatrixLit[",vecpos(3*parti+2,3*partj+2),"] * ",strings_o(modi-3),"*",strings_o(modj-3)," + ");
                fprintf(fid1,'%s%i%s%s%s%s%s',"MatrixLit[",vecpos(3*parti+2,3*partj+3),"] * ",strings_o(modi-3),"*",strings_o(modj)," + ");
                fprintf(fid1,'%s%i%s%s%s%s%s',"MatrixLit[",vecpos(3*parti+3,3*partj+2),"] * ",strings_o(modi),"*",strings_o(modj-3)," + ");
                fprintf(fid1,'%s%i%s%s%s%s%s\n',"MatrixLit[",vecpos(3*parti+3,3*partj+3),"] * ",strings_o(modi),"*",strings_o(modj),";");
            end
        end
    end
end