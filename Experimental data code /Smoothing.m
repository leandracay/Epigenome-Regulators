function[smoothedmCherry,PR2_smoothed,curvature2_PR2smoothed]=Smoothing(mCherrySubBackground)


PR=gradient(mCherrySubBackground,20);
[mCherrySubBackground,PR]=DerFilter(mCherrySubBackground,PR);

%same sized window
smoothedmCherry=smoothdata(mCherrySubBackground,2,'movmedian',10,'omitnan');
PR2=gradient(smoothedmCherry,20);
PR2_smoothed=smoothdata(PR2,2,'movmedian',10,'omitnan');
PR2_smoothed(isnan(mCherrySubBackground))=NaN;

curvature=gradient(PR2,20);
curvature2_PR2smoothed=gradient(PR2_smoothed,20);
curvature2_PR2smoothed(isnan(PR2_smoothed))=NaN;


end
