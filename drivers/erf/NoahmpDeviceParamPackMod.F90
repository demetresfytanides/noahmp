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

    nsoil = noahmp%config%domain%NumSoilLayer

    ! ---------------------------------------------------------------------
    ! PACK ORDER CONTRACT (mirror in ERF_NOAHMP_DeviceState.H unpack)
    ! [slice] scalars then per-layer soil arrays.
    ! ---------------------------------------------------------------------
    p = 0

    ! --- domain / config scalars ---
    p=p+1; buf(p) = real(noahmp%config%domain%NumSoilLayer,     kind_noahmp) ! [1] nsoil
    p=p+1; buf(p) = real(noahmp%config%domain%VegType,          kind_noahmp) ! [2] veg_type
    p=p+1; buf(p) = real(noahmp%config%domain%NumSwRadBand,     kind_noahmp) ! [3] num_sw_rad_band

    ! --- scalar energy param (table) ---
    p=p+1; buf(p) = noahmp%energy%param%SoilHeatCapacity                     ! [4] CSOIL_TABLE

    ! --- veg-indexed energy param ---
    p=p+1; buf(p) = noahmp%energy%param%HeightCanopyTop                      ! [5] HVT_TABLE(VegType)

    ! --- per-layer domain soil-depth array (nsoil) ---
    do L = 1, nsoil
       p=p+1; buf(p) = noahmp%config%domain%DepthSoilLayer(L)               ! [6..] ZSOIL
    enddo

    ! --- per-layer soil-type-indexed water param (nsoil) ---
    do L = 1, nsoil
       p=p+1; buf(p) = noahmp%water%param%SoilMoistureSat(L)                ! [..] SMCMAX_TABLE(SoilType(L))
    enddo

    nused = p

  end subroutine NoahmpDeviceParamPackColumn

end module NoahmpDeviceParamPackMod
