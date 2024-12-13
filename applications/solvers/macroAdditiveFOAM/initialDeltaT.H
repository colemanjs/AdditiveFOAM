scalar endTime = runTime.endTime().value();

forAll(sources.sources_, i)
{
    endTime = min(endTime, sources.sources_[i].beam().endTime());
}

scalar nPathIntervals =
    runTime.controlDict().lookupOrDefault<scalar>("nPathIntervals", 10.0);

const scalar deltaT0 =
    (endTime - runTime.startTime().value()) / nPathIntervals;

runTime.setDeltaT(deltaT0);
