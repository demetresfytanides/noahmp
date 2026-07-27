module NoahmpDeviceParamPackMod

!!! Pack the STATIC per-column Noah-MP parameters/config for the device (Kokkos GPU)
!!! land driver into a flat C-interop buffer, using the SAME validated *VarInTransfer
!!! lookups the CPU driver uses -- so the device params are correct by construction.
!!!
!!! Motivation: the ~130 device-state params the Kokkos NoahmpLandMain consumes are
!!! resolved in Fortran from ~300 module-level *_TABLE arrays (loaded from
!!! NoahmpTable.TBL) indexed by veg/soil/crop type. Those tables are NOT in the
!!! C++ NoahmpIO ABI mirror, so C++ cannot re-derive the params. Instead we run the
!!! transfer chain here (per column, once, since params are static for the run) and
!!! serialize the resolved params into a buffer C++ unpacks into the device LandState.
!!!
!!! The pack ORDER is a fixed contract shared with the C++ unpack in
!!! ERF_NOAHMP_DeviceState.H (erf_noahmp_gpu::unpack_column_params). Any change to
!!! the order/count on one side MUST be mirrored on the other; the nused return +
!!! the C++ NOAHMP_DEV_PARAM_COUNT guard catch drift at run time.
!!!
!!! Phase 2b VERTICAL SLICE: a representative subset spanning every param category
!!! (domain int, scalar table param, veg-indexed param, soil-indexed per-layer
!!! param, domain soil-layer array) to lock the Fortran->C++->LandState plumbing.
!!! Full ~130-field coverage is filled in incrementally after the slice builds green.

  use Machine
  use NoahmpIOVarType
  use NoahmpVarType
  use ConfigVarInitMod
  use ConfigVarInTransferMod
  use ForcingVarInitMod
  use ForcingVarInTransferMod
  use EnergyVarInitMod
  use EnergyVarInTransferMod
  use WaterVarInitMod
  use WaterVarInTransferMod
  use BiochemVarInitMod
  use BiochemVarInTransferMod

  implicit none

contains

  ! Resolve the static params for column (I,J) and serialize them into buf.
  !   buf(1:nbuf) : output buffer (kind_noahmp), caller-sized
  !   nused       : number of slots actually written (<= nbuf); C++ asserts == count
  ! Assumes the driver prelude scalars the transfers depend on are already set on
  ! NoahmpIO (ZSOIL, YEARLEN); NoahmpDeviceParamPrep sets them once per box.
  subroutine NoahmpDeviceParamPackColumn(noahmp, NoahmpIO, I, J, buf, nbuf, nused)

    implicit none

    type(noahmp_type),    intent(inout) :: noahmp
    type(NoahmpIO_type),  intent(inout) :: NoahmpIO
    integer,              intent(in)    :: I, J
    integer,              intent(in)    :: nbuf
    real(kind=kind_noahmp), intent(out) :: buf(nbuf)
    integer,              intent(out)   :: nused

    integer :: p          ! running write cursor
    integer :: L, nsoil

    NoahmpIO%I = I
    NoahmpIO%J = J

    ! Run the validated transfer chain into the 1-D column struct (same order as
    ! NoahmpDriverMain; params come out identical to the CPU path).
    call ConfigVarInitDefault  (noahmp)
    call ConfigVarInTransfer   (noahmp, NoahmpIO)
    call ForcingVarInitDefault (noahmp)
    call ForcingVarInTransfer  (noahmp, NoahmpIO)
    call EnergyVarInitDefault  (noahmp)
    call EnergyVarInTransfer   (noahmp, NoahmpIO)
    call WaterVarInitDefault   (noahmp)
    call WaterVarInTransfer    (noahmp, NoahmpIO)
    call BiochemVarInitDefault (noahmp)
    call BiochemVarInTransfer  (noahmp, NoahmpIO)

    nsoil = noahmp%config%domain%NumSoilLayer

    ! =====================================================================
    ! PACK ORDER CONTRACT (mirror EXACTLY in ERF_NOAHMP_DeviceState.H unpack).
    ! Grouped: config/domain -> energy scalar -> energy veg-indexed ->
    ! energy banded(nband) -> monthly LAI/SAI(12) -> water scalar ->
    ! water veg/slope-indexed -> per-layer soil arrays(nsoil).
    ! Count = noahmp_dev_param_count(nsoil,nband) on both sides.
    ! =====================================================================
    p = 0
    associate(dm => noahmp%config%domain, nml => noahmp%config%nmlist, &
              ep => noahmp%energy%param,   wp  => noahmp%water%param)

    ! ---- [G1] config / domain scalars (ints stored as reals) ----
    p=p+1; buf(p) = real(dm%NumSoilLayer,   kind_noahmp)
    p=p+1; buf(p) = real(dm%NumSwRadBand,   kind_noahmp)
    p=p+1; buf(p) = real(dm%VegType,        kind_noahmp)
    p=p+1; buf(p) = real(dm%SurfaceType,    kind_noahmp)
    p=p+1; buf(p) = real(dm%IndexIcePoint,  kind_noahmp)
    p=p+1; buf(p) = real(dm%IndexBarrenPoint,kind_noahmp)
    p=p+1; buf(p) = real(dm%IndexWaterPoint,kind_noahmp)
    p=p+1; buf(p) = real(dm%CropType,       kind_noahmp)
    p=p+1; buf(p) = real(merge(1,0,dm%FlagUrban), kind_noahmp)
    p=p+1; buf(p) = real(wp%NumSoilLayerRoot, kind_noahmp)
    ! namelist options
    p=p+1; buf(p) = real(nml%OptDynamicVeg,            kind_noahmp)
    p=p+1; buf(p) = real(nml%OptRainSnowPartition,     kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSoilWaterTranspiration,kind_noahmp)
    p=p+1; buf(p) = real(nml%OptGroundResistanceEvap,  kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSurfaceDrag,           kind_noahmp)
    p=p+1; buf(p) = real(nml%OptStomataResistance,     kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSnowAlbedo,            kind_noahmp)
    p=p+1; buf(p) = real(nml%OptCanopyRadiationTransfer,kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSnowSoilTempTime,      kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSnowThermConduct,      kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSoilTemperatureBottom, kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSoilSupercoolWater,    kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSoilPermeabilityFrozen,kind_noahmp)
    p=p+1; buf(p) = real(nml%OptTileDrainage,          kind_noahmp)
    p=p+1; buf(p) = real(nml%OptRunoffSurface,         kind_noahmp)
    p=p+1; buf(p) = real(nml%OptRunoffSubsurface,      kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSnowCompaction,        kind_noahmp)
    p=p+1; buf(p) = real(nml%OptWetlandModel,          kind_noahmp)
    p=p+1; buf(p) = real(nml%OptSnowCoverGround,       kind_noahmp)
    p=p+1; buf(p) = real(nml%OptCropModel,             kind_noahmp)
    ! domain reals
    p=p+1; buf(p) = dm%MainTimeStep
    p=p+1; buf(p) = dm%SoilTimeStep
    p=p+1; buf(p) = real(dm%NumSoilTimeStep, kind_noahmp)
    p=p+1; buf(p) = dm%GridSize
    p=p+1; buf(p) = dm%DepthSoilTempBottom
    p=p+1; buf(p) = dm%CosSolarZenithAngle

    ! ---- [G2] energy scalar params ----
    p=p+1; buf(p) = ep%SoilHeatCapacity
    p=p+1; buf(p) = ep%SnowAgeFacBats
    p=p+1; buf(p) = ep%SnowGrowVapFacBats
    p=p+1; buf(p) = ep%SnowSootFacBats
    p=p+1; buf(p) = ep%SnowGrowFrzFacBats
    p=p+1; buf(p) = ep%SolarZenithAdjBats
    p=p+1; buf(p) = ep%FreshSnoAlbVisBats
    p=p+1; buf(p) = ep%FreshSnoAlbNirBats
    p=p+1; buf(p) = ep%SnoAgeFacDifVisBats
    p=p+1; buf(p) = ep%SnoAgeFacDifNirBats
    p=p+1; buf(p) = ep%SzaFacDirVisBats
    p=p+1; buf(p) = ep%SzaFacDirNirBats
    p=p+1; buf(p) = ep%UpscatterCoeffSnowDir
    p=p+1; buf(p) = ep%UpscatterCoeffSnowDif
    p=p+1; buf(p) = ep%EmissivitySnow
    p=p+1; buf(p) = ep%EmissivitySoilLake(1)   ! (:) idx 1=soil,2=lake -> soil path
    p=p+1; buf(p) = ep%EmissivityIceSfc
    p=p+1; buf(p) = ep%RoughLenMomSnow
    p=p+1; buf(p) = ep%RoughLenMomSoil
    p=p+1; buf(p) = ep%RoughLenMomLake
    p=p+1; buf(p) = ep%ResistanceSoilExp
    p=p+1; buf(p) = ep%ResistanceSnowSfc
    p=p+1; buf(p) = ep%VegFracAnnMax
    p=p+1; buf(p) = ep%VegFracGreen
    p=p+1; buf(p) = wp%SnowMassFullCoverOld   ! SWEMX_TABLE (Water param; snow-aging in)

    ! ---- [G3] energy veg-indexed params ----
    p=p+1; buf(p) = ep%TreeCrownRadius
    p=p+1; buf(p) = ep%HeightCanopyTop
    p=p+1; buf(p) = ep%HeightCanopyBot
    p=p+1; buf(p) = ep%RoughLenMomVeg
    p=p+1; buf(p) = ep%CanopyWindExtFac
    p=p+1; buf(p) = ep%TreeDensity
    p=p+1; buf(p) = ep%CanopyOrientIndex
    p=p+1; buf(p) = ep%HeatCapacCanFac

    ! ---- [G4] energy banded params (1:NumSwRadBand) ----
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%ReflectanceLeaf(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%ReflectanceStem(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%TransmittanceLeaf(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%TransmittanceStem(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%AlbedoSoilSat(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%AlbedoSoilDry(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%AlbedoLakeFrz(L)
    enddo
    do L = 1, dm%NumSwRadBand
       p=p+1; buf(p) = ep%ScatterCoeffSnow(L)
    enddo

    ! ---- [G5] monthly LAI/SAI (1:12) ----
    do L = 1, 12
       p=p+1; buf(p) = ep%LeafAreaIndexMon(L)
    enddo
    do L = 1, 12
       p=p+1; buf(p) = ep%StemAreaIndexMon(L)
    enddo

    ! ---- [G6] water scalar params ----
    p=p+1; buf(p) = wp%SnowCompactBurdenFac
    p=p+1; buf(p) = wp%SnowCompactAgingFac1
    p=p+1; buf(p) = wp%SnowCompactAgingFac2
    p=p+1; buf(p) = wp%SnowCompactAgingFac3
    p=p+1; buf(p) = wp%SnowCompactAgingMax
    p=p+1; buf(p) = wp%SnowViscosityCoeff
    p=p+1; buf(p) = wp%SnowLiqFracMax
    p=p+1; buf(p) = wp%SnowLiqHoldCap
    p=p+1; buf(p) = wp%SnowLiqReleaseFac
    p=p+1; buf(p) = wp%SoilConductivityRef
    p=p+1; buf(p) = wp%SoilInfilFacRef
    p=p+1; buf(p) = wp%GroundFrzCoeff
    p=p+1; buf(p) = wp%SnowfallDensityMax
    p=p+1; buf(p) = wp%SnowMassFullCoverOld
    p=p+1; buf(p) = wp%SoilMatPotentialWilt
    p=p+1; buf(p) = wp%SnoWatEqvMaxGlacier
    p=p+1; buf(p) = wp%SoilInfilMaxCoeff
    p=p+1; buf(p) = wp%SoilImpervFracCoeff
    p=p+1; buf(p) = wp%SoilDrainSlope

    ! ---- [G7] water veg-indexed params ----
    p=p+1; buf(p) = wp%CanopyLiqHoldCap
    p=p+1; buf(p) = wp%SnowMeltFac
    p=p+1; buf(p) = wp%SnowCoverFac

    ! ---- [G8] per-layer soil arrays (1:nsoil) ----
    do L = 1, nsoil
       p=p+1; buf(p) = dm%DepthSoilLayer(L)
    enddo
    ! (ThicknessSoilLayer is NOT packed: unresolved by the transfer chain
    !  until GeneralInit; the device general_init derives it from depth.)
    do L = 1, nsoil
       p=p+1; buf(p) = ep%SoilQuartzFrac(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilMoistureSat(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilMoistureWilt(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilMoistureFieldCap(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilMoistureDry(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilWatDiffusivitySat(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilWatConductivitySat(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilExpCoeffB(L)
    enddo
    do L = 1, nsoil
       p=p+1; buf(p) = wp%SoilMatPotentialSat(L)
    enddo

    ! ---- [G9] vegetated flux tile: Ball-Berry stomata + Jarvis + transpiration
    !      params (energy + biochem). Consumed only when the veg tile activates. ----
    associate(bp => noahmp%biochem%param)
    p=p+1; buf(p) = ep%ConductanceLeafMin        ! BP_TABLE
    p=p+1; buf(p) = ep%Co2MmConst25C             ! KC25_TABLE
    p=p+1; buf(p) = ep%O2MmConst25C              ! KO25_TABLE
    p=p+1; buf(p) = ep%Co2MmConstQ10             ! AKC_TABLE
    p=p+1; buf(p) = ep%O2MmConstQ10              ! AKO_TABLE
    p=p+1; buf(p) = ep%RadiationStressFac        ! RGL_TABLE (Jarvis)
    p=p+1; buf(p) = ep%ResistanceStomataMin      ! RS_TABLE
    p=p+1; buf(p) = ep%ResistanceStomataMax      ! RSMAX_TABLE
    p=p+1; buf(p) = ep%AirTempOptimTransp        ! TOPT_TABLE
    p=p+1; buf(p) = ep%VaporPresDeficitFac       ! HS_TABLE
    p=p+1; buf(p) = ep%LeafDimLength             ! DLEAF_TABLE (r_leaf)
    p=p+1; buf(p) = ep%CanopyWindExtFac          ! CWPVT_TABLE (r_leaf)
    p=p+1; buf(p) = bp%QuantumEfficiency25C      ! QE25_TABLE
    p=p+1; buf(p) = bp%CarboxylRateMax25C        ! VCMX25_TABLE
    p=p+1; buf(p) = bp%CarboxylRateMaxQ10        ! AVCMX_TABLE
    p=p+1; buf(p) = bp%PhotosynPathC3            ! C3PSN_TABLE
    p=p+1; buf(p) = bp%SlopeConductToPhotosyn    ! MP_TABLE
    p=p+1; buf(p) = bp%NitrogenConcFoliageMax    ! FOLNMX_TABLE
    p=p+1; buf(p) = noahmp%energy%state%PressureAtmosCO2  ! CO2_TABLE*PressAir
    p=p+1; buf(p) = noahmp%energy%state%PressureAtmosO2   ! O2_TABLE *PressAir
    end associate

    end associate
    nused = p

  end subroutine NoahmpDeviceParamPackColumn

end module NoahmpDeviceParamPackMod
